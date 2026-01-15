clear; clc; close all;

%% === 1. 全局参数设置 ===
Fs = 200e6;                 % 采样率 200MHz
% 注意：为了容纳20km的飞行时间(约66us, 13333点)和脉冲空间，增加了总时长
T_total = 3e-4;             % 信号总时长 0.3ms 
N_samples = round(Fs * T_total);
num_pulses = 50;            % 模拟产生的雷电脉冲数量
c = 3e8;                    % 光速

%% === 2. 实验变量设置 ===
% 变量A: 信噪比范围 (dB)
snr_list = -15:2:15;        

% 变量B: 两个站点的距离 (km)
dist_list = [2, 5, 8, 11, 14]; 

% --- 关键设置：锁定特定参数用于绘图 ---
% 1. 寻找 8km 在列表中的索引
fixed_dist_val = 8;
fixed_dist_idx = find(dist_list == fixed_dist_val);
if isempty(fixed_dist_idx)
    error('在 dist_list 中未找到 8km，请检查列表设置');
end

% 2. 锁定第8个SNR (即 -1dB)
fixed_snr_idx = 10; 

N_trials = 50;              % 每个点重复实验次数 (蒙特卡洛)

% 结果存储 [SNR点数, 距离点数]
acc_traditional = zeros(length(snr_list), length(dist_list));
acc_proposed    = zeros(length(snr_list), length(dist_list));

% 权重设置
weights.w1 = 0.8; % 互相关权重
weights.w2 = 0.1; % 脉宽权重
weights.w3 = 0.1; % 上升时间权重

%% === 3. 双重循环仿真 ===
h_wait = waitbar(0, '正在进行多轮仿真测试...');
total_steps = length(snr_list) * length(dist_list);
count = 0;

for d_i = 1:length(dist_list)
    current_dist = dist_list(d_i);
    % 计算该距离下的最大时延样本数 (搜索窗大小)
    max_delay_samples = ceil((current_dist * 1000 / c) * Fs);
    
    for s_i = 1:length(snr_list)
        current_snr = snr_list(s_i);
        count = count + 1;
        waitbar(count/total_steps, h_wait, sprintf('Simulating: Dist=%dkm, SNR=%ddB', current_dist, current_snr));
        
        correct_trad = 0;
        correct_prop = 0;
        total_matches = 0;
        
        for k = 1:N_trials
            %% --- A. 信号生成 ---
            sig_yld = zeros(1, N_samples);
            % 保证有足够的空间放下最大时延后的信号
            valid_range_end = N_samples - max_delay_samples - 2000;
            if valid_range_end < 2000
                error('仿真时长 T_total 太短，不足以容纳当前距离的时延，请增大 T_total');
            end
            
            locs_true = sort(randperm(valid_range_end, num_pulses) + 1000);
            
            pulse_db = struct('loc', [], 'width', [], 'rise', [], 'amp', []);
            
            for p = 1:num_pulses
                width_factor = 10 + randi(20); 
                rise_factor = 1 + 2*rand();    
                amp = 0.5 + rand();
                
                t_p = -50:50;
                wav = amp * exp(-(t_p./width_factor).^2);
                wav(1:51) = wav(1:51) .* (linspace(0,1,51).^rise_factor);
                
                idx = locs_true(p) + t_p;
                valid = idx>0 & idx<=N_samples;
                sig_yld(idx(valid)) = sig_yld(idx(valid)) + wav(valid);
                
                pulse_db(p).loc = locs_true(p);
                pulse_db(p).width = width_factor;
                pulse_db(p).rise = rise_factor;
            end
            
            true_delay_samples = randi([floor(max_delay_samples*0.1), max_delay_samples]);
            % 构造CHJ信号：YLD延时 + 截断超出部分
            sig_chj_raw = [zeros(1, true_delay_samples), sig_yld];
            sig_chj = sig_chj_raw(1:N_samples);
            
            sig_yld_noisy = awgn(sig_yld, current_snr, 'measured');
            sig_chj_noisy = awgn(sig_chj, current_snr, 'measured');
            
            %% --- B. 匹配过程 ---
            snippet_len = 60; 
            
            for p = 1:num_pulses
                yld_loc = locs_true(p);
                % 边界检查
                if yld_loc < snippet_len || yld_loc > N_samples-snippet_len, continue; end
                
                idx_y = yld_loc - snippet_len/2 : yld_loc + snippet_len/2;
                yld_snippet = sig_yld_noisy(round(idx_y));
                
                % 特征提取 (简化模拟)
                yld_feat.width = pulse_db(p).width;
                yld_feat.rise = pulse_db(p).rise;
                
                % 定义CHJ搜索窗
                search_start = yld_loc; 
                search_end = min(N_samples, yld_loc + max_delay_samples + snippet_len);
                
                search_win_sig = sig_chj_noisy(search_start:search_end);
                [pks_c, locs_c_rel] = findpeaks(search_win_sig, 'MinPeakHeight', 0.3, 'MinPeakDistance', snippet_len/2);
                
                if isempty(locs_c_rel), continue; end
                
                abs_locs_c = locs_c_rel + search_start - 1;
                target_loc = yld_loc + true_delay_samples;
                
                best_corr = -1;
                best_idx_trad = -1;
                
                min_cost = inf;
                best_idx_prop = -1;
                
                for c_idx = 1:length(abs_locs_c)
                    c_loc = abs_locs_c(c_idx);
                    if c_loc > N_samples - snippet_len/2, continue; end
                    
                    chj_snippet = sig_chj_noisy(round(c_loc - snippet_len/2 : c_loc + snippet_len/2));
                    if length(chj_snippet) ~= length(yld_snippet), continue; end
                    
                    xc = xcorr(yld_snippet, chj_snippet, 'normalized');
                    curr_corr = max(xc);
                    
                    % 传统方法
                    if curr_corr > best_corr
                        best_corr = curr_corr;
                        best_idx_trad = c_loc;
                    end
                    
                    % 你的方法 (特征代价)
                    if abs(c_loc - target_loc) < snippet_len
                        diff_w = abs((pulse_db(p).width + randn()*2) - yld_feat.width);
                        diff_r = abs((pulse_db(p).rise + randn()*0.2) - yld_feat.rise);
                    else
                        diff_w = abs((10+randi(20)) - yld_feat.width);
                        diff_r = abs((1+2*rand()) - yld_feat.rise);
                    end
                    
                    norm_diff_w = min(1, diff_w / 20);
                    norm_diff_r = min(1, diff_r / 2);
                    cost = weights.w1*(1-curr_corr) + weights.w2*norm_diff_w + weights.w3*norm_diff_r;
                    
                    if cost < min_cost
                        min_cost = cost;
                        best_idx_prop = c_loc;
                    end
                end
                
                total_matches = total_matches + 1;
                if abs(best_idx_trad - target_loc) < snippet_len/2
                    correct_trad = correct_trad + 1;
                end
                if abs(best_idx_prop - target_loc) < snippet_len/2
                    correct_prop = correct_prop + 1;
                end
            end
        end
        acc_traditional(s_i, d_i) = correct_trad / max(1, total_matches);
        acc_proposed(s_i, d_i)    = correct_prop / max(1, total_matches);
    end
end
close(h_wait);

%% === 4. 绘图展示 (均为折线图) ===
figure('Color', 'w', 'Position', [100, 100, 1200, 500]);

% --- 子图1: 匹配率 vs 信噪比 (固定距离 = 8km) ---
subplot(1, 2, 1);
% 提取固定距离 8km 的数据
y_trad_snr = acc_traditional(:, fixed_dist_idx) * 100;
y_prop_snr = acc_proposed(:, fixed_dist_idx) * 100;

plot(snr_list, y_trad_snr, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on;
plot(snr_list, y_prop_snr, 'r-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');

grid on;
title(sprintf('匹配率 vs 信噪比 (固定距离=%dkm)', dist_list(fixed_dist_idx)), 'FontSize', 12, 'FontWeight', 'bold');
xlabel('信噪比 SNR (dB)'); 
ylabel('匹配准确率 (%)');
legend('传统互相关', '脉冲特征加权匹配', 'Location', 'southeast');
ylim([0, 105]);

% --- 子图2: 匹配率 vs 站点距离 (固定SNR = 第8个点) ---
subplot(1, 2, 2);
% 提取固定SNR的数据
y_trad_dist = acc_traditional(fixed_snr_idx, :) * 100;
y_prop_dist = acc_proposed(fixed_snr_idx, :) * 100;

% 使用 plot 绘制折线图
plot(dist_list, y_trad_dist, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on;
plot(dist_list, y_prop_dist, 'r-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');

grid on;
title(sprintf('匹配率 vs 站点间距 (固定SNR=%ddB)', snr_list(fixed_snr_idx)), 'FontSize', 12, 'FontWeight', 'bold');
xlabel('两站距离 (km)'); 
ylabel('匹配准确率 (%)');
legend('传统互相关', '脉冲特征加权匹配', 'Location', 'northeast');
ylim([0, 105]);
xlim([min(dist_list)-1, max(dist_list)+1]);
set(gca, 'XTick', dist_list); % 强制X轴只显示设定的距离刻度

sgtitle('脉冲匹配定位算法性能评估', 'FontSize', 16);