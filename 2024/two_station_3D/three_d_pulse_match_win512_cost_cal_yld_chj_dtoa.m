%% ========================================================================
%  主程序: 脉冲级特征匹配定位算法
%  ========================================================================

clear; clc; close all;

%% --- 0. 初始化和参数定义 ---
yld_sit = [0, 0, 0];
chj_sit = [2003.7972, -7844.7836, -27];
% yld相对于chj的位置
p = chj_sit-yld_sit;
dist = 8.0967e3; %单位：米
c = 0.299792458;
W = 5e8; % 时间误差
match_signal_length = 6000;
% 两个站点的距离除以光速
rough_dtoa_samples = 6000;
step = 1e6;
signal_length=step;
start_signal_loc = 3.6e8;
mapping_start_signal_loc = 3.6e8;
end_signal_loc = 5.6e8;

% 所有信号的开始位置
all_start_signal_loc = start_signal_loc:step:end_signal_loc;
all_match_results = []; % 存储DTOA优化后的匹配信息
% YLD站内天线局部坐标
yld_ant_local_coords = [
    0,   0,   0;  % 天线1
    -23.3246,  -8.8824,   0;  % 天线2
    -31.7222,  14.6305,   0   % 天线3
    ];
% CHJ站内天线局部坐标
chj_ant_local_coords = [
    0,   0,   0;  % 天线1
    -2.0622, 41.5985,   0;  % 天线2
    28.4316, 23.5237,   0   % 天线3
    ];

% 将局部天线坐标转换为全局坐标
yld_ant_global_coords = yld_ant_local_coords + yld_sit;
chj_ant_global_coords = chj_ant_local_coords + chj_sit;

baseline_ant_indices = [
    1, 2; % 基线12
    1, 3; % 基线13
    2, 3  % 基线23
    ];
% dtoa优化时的搜索范围，只会在初始点的上下左右10范围内进行优化，需要调整
search_radius_xy = 100;
search_radius_z  = 100;
% --- 所有基线的天线全局坐标 (6条站内 + 1条站间) ---
all_baseline_definitions = zeros(7, 6);
all_baseline_definitions(1, :) = [yld_ant_global_coords(1,:), yld_ant_global_coords(2,:)];
all_baseline_definitions(2, :) = [yld_ant_global_coords(1,:), yld_ant_global_coords(3,:)];
all_baseline_definitions(3, :) = [yld_ant_global_coords(2,:), yld_ant_global_coords(3,:)];
all_baseline_definitions(4, :) = [chj_ant_global_coords(1,:), chj_ant_global_coords(2,:)];
all_baseline_definitions(5, :) = [chj_ant_global_coords(1,:), chj_ant_global_coords(3,:)];
all_baseline_definitions(6, :) = [chj_ant_global_coords(2,:), chj_ant_global_coords(3,:)];
% 站间基线 (YLD参考天线 - CHJ参考天线)
all_baseline_definitions(7, :) = [yld_sit, chj_sit];
% 记录所有优化后的残差 = 理论DTOA - 测量DTOA
all_residuals_history = [];
% --- DTOA测量不确定度  ---  站间基线的不确定度需要尝试，如果设置成5000几乎与初始点一致
sigmas_ns = [
    2; 2; 2; % YLD 站内基线 DTOA 不确定度 (ns)
    2; 2; 2; % CHJ 站内基线 DTOA 不确定度 (ns)
    500         % 站间基线 DTOA 不确定度 (ns)
    ];



% 信号和算法参数
fs = 200e6;
ts_ns = 1 / fs * 1e9;

noise_analysis_length = 1e5;
threshold_std_multiplier = 4;
noise_chj = read_signal('..\\2024 822 85933.651462CH1.dat', noise_analysis_length, noise_analysis_length);
filtered_noise_chj = filter_bp(noise_chj, 30e6, 80e6, 5);
threshold_chj = mean(filtered_noise_chj) + threshold_std_multiplier * std(filtered_noise_chj);
noise_yld = read_signal('..\\20240822165932.6610CH1.dat', noise_analysis_length, noise_analysis_length);
filtered_noise_yld = filter_bp(noise_yld, 30e6, 80e6, 5);
threshold_yld = mean(filtered_noise_yld) + threshold_std_multiplier * std(filtered_noise_yld);

snippet_len = 512;     % 脉冲片段长度
% 代价函数权重
weights.w1 = 0.8; % 波形相似度权重
weights.w2 = 0.1; % FWHM 权重 脉冲宽度
weights.w3 = 0.1; % 上升时间权重
cost_threshold = 5; % 只有总代价小于此阈值的匹配才被接受, 需要调整

% DTOA时间搜索窗
time_window_samples = 6000;

for j = 1:numel(all_start_signal_loc)-1
    start_read_loc_yld = all_start_signal_loc(j);
    start_read_loc_chj = start_read_loc_yld;
    end_read_loc_yld = all_start_signal_loc(j+1);
    if  start_read_loc_yld < mapping_start_signal_loc || end_read_loc_yld > end_signal_loc
        continue
    end
    % 记录处理的位置
    fprintf('正在处理的信号位置：%d -- %d \n', start_read_loc_yld, end_read_loc_yld);
    yld_ch1 =read_signal('..\\20240822165932.6610CH1.dat',signal_length,start_read_loc_yld);
    yld_ch2 =read_signal('..\\20240822165932.6610CH2.dat',signal_length,start_read_loc_yld);
    yld_ch3 =read_signal('..\\20240822165932.6610CH3.dat',signal_length,start_read_loc_yld);
    filtered_yld_signal1 = filter_bp(yld_ch1,30e6,80e6,5);
    filtered_yld_signal2 = filter_bp(yld_ch2,30e6,80e6,5);
    filtered_yld_signal3 = filter_bp(yld_ch3,30e6,80e6,5);
    %% --- 1. 独立脉冲编目  ---
    % 创建脉冲目录
    yld_catalog = create_pulse_catalog(filtered_yld_signal1, filtered_yld_signal2, filtered_yld_signal3, fs, threshold_yld, snippet_len);
    % 检查目录是否为空
    if isempty(yld_catalog)
        fprintf('yld站的脉冲目录为空，无法继续。');
        continue
    end
    num_yld_pulses = numel(yld_catalog);
    h = waitbar(0, 'Processing...');
    for i =1 :num_yld_pulses
        waitbar(i/num_yld_pulses, h, sprintf('位置：%d -- %d ；进度： %.2f%%', start_read_loc_yld, end_read_loc_yld, i/num_yld_pulses*100));
        current_yld_pulse = yld_catalog(i);
        current_yld_pulse.loc = start_read_loc_yld - start_signal_loc + current_yld_pulse.loc;
        current_chj_pulse_loc = floor(current_yld_pulse.loc *5/5.00390073);
        chj_ch1 =read_signal('..\\2024 822 85933.651462CH1.dat',match_signal_length*2,start_signal_loc+ 34236596 + current_chj_pulse_loc - match_signal_length);
        chj_ch2 =read_signal('..\\2024 822 85933.651462CH2.dat',match_signal_length*2,start_signal_loc+ 34236596 + current_chj_pulse_loc - match_signal_length);
        chj_ch3 =read_signal('..\\2024 822 85933.651462CH3.dat',match_signal_length*2,start_signal_loc+ 34236596 + current_chj_pulse_loc - match_signal_length + 215/5);
        filtered_chj_signal1 = filter_bp(chj_ch1,30e6,80e6,5);
        filtered_chj_signal2 = filter_bp(chj_ch2,30e6,80e6,5);
        filtered_chj_signal3 = filter_bp(chj_ch3,30e6,80e6,5);
        %% Step2 根据引雷点的信号窗口得到匹配到的从化局的信号
        chj_catalog = create_pulse_catalog(filtered_chj_signal1, filtered_chj_signal2, filtered_chj_signal3, fs, threshold_chj, snippet_len);
        [matched_chj_pulse, match_cost] = find_best_match(current_yld_pulse, chj_catalog, weights);
        % --- 2.3 如果找到一个好的匹配，则进行精处理 ---
        if ~isempty(matched_chj_pulse) && match_cost < cost_threshold
            matched_chj_pulse.loc = matched_chj_pulse.loc + current_chj_pulse_loc - match_signal_length;
        else
            continue
        end
        % a. 获取YLD和CHJ的2D测向结果
        [yld_t12, yld_t13,yld_t23,yld_az, yld_el, yld_Rcorr, yld_t123] = get_2d_from_pulse(current_yld_pulse, 'yld');
        [chj_t12, chj_t13, chj_t23, chj_az, chj_el, chj_Rcorr, chj_t123] = get_2d_from_pulse(matched_chj_pulse, 'chj');
        chj_ch1 = matched_chj_pulse.snippet;
        yld_ch1 = current_yld_pulse.snippet;

        if length(chj_ch1)~= length(yld_ch1)
            continue
        end

        [r_gcc, lags_gcc] = xcorr(chj_ch1, yld_ch1, 'normalized');
        max_Rgcc = max(r_gcc);
        if yld_Rcorr < 0.6 || abs(yld_t123) > 1 || yld_el  > 80 || yld_el < 10 || chj_Rcorr < 0.3 || abs(chj_t123) > 1 || max_Rgcc < 0.1
            continue
        end

        [R1_x, R1_y, R1_z] = sph2cart(deg2rad(90-yld_az), deg2rad(yld_el),1);
        [R2_x, R2_y, R2_z] = sph2cart(deg2rad(90-chj_az), deg2rad(chj_el),1);
        A1 = [R1_x, R1_y, R1_z];
        A2 = [R2_x, R2_y, R2_z];
        C = cross(A1, A2);
        if norm(C) < eps
            continue;  % 避免除以零
        end
        c_unit = C  / norm(C);  % 单位向量
        M = [A1(1), -A2(1), c_unit(1);
            A1(2), -A2(2), c_unit(2);
            A1(3), -A2(3), c_unit(3)];
        % 使用克莱姆法则求R1,R2,R3的标量
        detM = det(M);
        detR1 = det([p', M(:,2), M(:,3)]);
        detR2 = det([M(:,1), p', M(:,3)]);
        detR3 = det([M(:,1), M(:,2), p']);
        R1_value = detR1 / detM;
        R2_value = detR2 / detM;
        R3_value = detR3 / detM;
        R1 = R1_value * A1;
        R2 = R2_value * A2;
        R3 = R3_value/norm(C)* C;
        if R1_value <= R2_value
            % 使用第一个公式
            sub_S = R1 + (R1_value / R2_value)*(R1_value / (R1_value + R2_value)) * R3;
        else
            % 使用第二个公式
            sub_S = R2 - (R2_value / R1_value)* (R2_value / (R1_value + R2_value))  * R3 + p;
        end
        if sub_S(3) < 0
            sub_S = -sub_S;
        end
        if ~isempty(sub_S)
            t_chj = sqrt(sum((sub_S - chj_sit).^2))/c;
            t_yld = sqrt(sum((sub_S - yld_sit).^2))/c;
            dlta_t = t_yld-t_chj;
            dlta_T = (current_yld_pulse.loc-matched_chj_pulse.loc/5*5.00390073)*5;
            dlta = abs(dlta_t-dlta_T);
            if dlta <= W
                %% Step 5: 差分到达时间 (DTOA) 技术
                % --- 收集所有测量的DTOA值 (单位: ns) ---

                all_measured_dtoas_ns = zeros(7,1);
                all_measured_dtoas_ns(1:3) = [yld_t12,yld_t13,yld_t23];  % YLD站内 (ns)
                all_measured_dtoas_ns(4:6) = [chj_t12,chj_t13,chj_t23];    % CHJ站内 (ns)
                % 站间DTOA (ns)
                all_measured_dtoas_ns(7) = dlta_T;
                % --- DTOA 优化 ---
                S_initial = sub_S;   % 使用三角测量结果作为初值
                objective_fun = @(S_vec) dtoa_objective_function_ns(S_vec, ...
                    all_baseline_definitions, ...
                    all_measured_dtoas_ns, ...
                    sigmas_ns, ...
                    c);


                dynamic_lb = S_initial - [search_radius_xy, search_radius_xy, search_radius_z];
                dynamic_ub = S_initial + [search_radius_xy, search_radius_xy, search_radius_z];
                options = optimoptions('lsqnonlin', 'TolFun', 1e-6, 'MaxIter', 100);
                S_optimized = S_initial; % 默认值
                optimization_successful = false;
                try
                    [S_optimized_temp, resnorm, ~, exitflag, ~] = lsqnonlin(objective_fun, S_initial, dynamic_lb, dynamic_ub, options);
                    if exitflag > 0 % 检查优化是否成功收敛
                        S_optimized = S_optimized_temp;
                        optimization_successful = true;
                    else
                        fprintf('优化未收敛 (exitflag=%d) for YLD event %d (sample loc %d)\n', exitflag, i, current_yld_start_abs_samples);
                    end
                catch ME_optim
                    fprintf('DTOA优化错误 for YLD event %d (sample loc %d): %s\n', i, yld_start_loc(i), ME_optim.message);
                end

                if optimization_successful
                    % 用最终的优化解 S_optimized 再次调用目标函数，获取最终的理论DTOA
                    final_theoretical_dtoas_ns = dtoa_objective_function_ns(S_optimized, ...
                        all_baseline_definitions, ...
                        all_measured_dtoas_ns, ...
                        sigmas_ns, ...
                        c);
                    % 残差 = 理论DTOA - 测量DTOA
                    final_residuals_ns = final_theoretical_dtoas_ns;
                    chi_square_red = sum(final_residuals_ns.^2) / (length(final_residuals_ns) - 3);
                    all_residuals_history = [all_residuals_history; final_residuals_ns(:)'];

                end

                % --- 使用优化后的S重新进行时间校验 ---
                % 理论传播时间 (ns)
                t_chj_optimized_ns = norm(S_optimized - chj_sit) / c;
                t_yld_optimized_ns = norm(S_optimized - yld_sit) / c;
                % 理论站间DTOA (ns),
                theoretical_inter_dtoa_optimized_ns = t_yld_optimized_ns - t_chj_optimized_ns;

                % 最终校验差 (ns)
                final_dlta_check_ns = abs(theoretical_inter_dtoa_optimized_ns - dlta_T);

                if final_dlta_check_ns <= W && optimization_successful
                    if S_optimized(3) < 0
                        S_optimized = -S_optimized;
                    end
                    %                     我们还可以记录损失
                    match_info_dtoa = struct(...
                        'yld_start_loc', start_signal_loc + current_yld_pulse.loc - snippet_len/2 + 1, ...
                        'chj_loc', matched_chj_pulse.loc + start_signal_loc + 34226222 - snippet_len/2 +1, ...
                        'chj_azimuth', chj_az, ...
                        'chj_elevation', chj_el, ...
                        'r_gccs', max_Rgcc, ...
                        'dlta', dlta, ...
                        'R3_value', R3_value, ...
                        'chi_square_red', chi_square_red,  ...
                        'cost', match_cost,  ...
                        'dlta_T',dlta_T, ...
                        'x_tri', S_initial(1), ...  % 【新增】明确保存三角测量 X
                        'y_tri', S_initial(2), ...  % 【新增】明确保存三角测量 Y
                        'z_tri', S_initial(3), ...  % 【新增】明确保存三角测量 Z
                        'x_dtoa', S_optimized(1), ... % 【修改】明确命名为 DTOA X
                        'y_dtoa', S_optimized(2), ... % 【修改】明确命名为 DTOA Y
                        'z_dtoa', S_optimized(3) ...  % 【修改】明确命名为 DTOA Z
                        );
                    all_match_results = [all_match_results; match_info_dtoa];
                else
                    % 如果DTOA优化失败或校验不通过，保存三角测量结果
                    if dlta <= W
                        if S_initial(3) < 0
                            S_initial = -S_initial;
                        end
                        match_info_dtoa = struct(...
                            'yld_start_loc', start_signal_loc + current_yld_pulse.loc - snippet_len/2 + 1, ...
                            'chj_loc', matched_chj_pulse.loc + start_signal_loc + 34226222 - snippet_len/2 +1, ...
                            'chj_azimuth', chj_az, ...
                            'chj_elevation', chj_el, ...
                            'r_gccs', max_Rgcc, ...
                            'dlta', dlta, ...
                            'R3_value', R3_value, ...
                            'chi_square_red', 1.5 , ...
                            'cost', match_cost,  ...
                            'dlta_T',dlta_T, ...
                            'x_tri', S_initial(1), ...  % 【新增】明确保存三角测量 X
                            'y_tri', S_initial(2), ...  % 【新增】明确保存三角测量 Y
                            'z_tri', S_initial(3), ...  % 【新增】明确保存三角测量 Z
                            'x_dtoa', S_initial(1), ... % 【修改】明确命名为 DTOA X
                            'y_dtoa', S_initial(2), ... % 【修改】明确命名为 DTOA Y
                            'z_dtoa', S_initial(3) ...  % 【修改】明确命名为 DTOA Z
                            );
                        all_match_results = [all_match_results; match_info_dtoa];
                    end
                end
            end
        end
    end
    close(h);
end




all_match_table = struct2table(all_match_results);

% 步骤2：用 writetable 保存为 CSV（支持 'Delimiter' 参数）
writetable(all_match_table, ...
           'results\3d_win512_cost_cal_yld_chj_dtoa_3.6e8_5.6e8_with_initial_results.csv', ...  % 文件名
           'Encoding', 'UTF-8', ...                    % 编码（确保中文正常）
           'Delimiter', ',', ...                       % 指定逗号分隔（CSV标准）
           'WriteVariableNames', true);                % 保存列名（结构体字段名）


% --- DTOA 目标函数 ---
function errors = dtoa_objective_function_ns(S_guess_xyz, all_baseline_ant_coords, all_measured_dtoas_ns, all_sigmas_ns, c_mpns)
num_baselines = size(all_baseline_ant_coords, 1);
errors = zeros(num_baselines, 1);
S_guess_xyz = S_guess_xyz(:)'; % 确保是行向量 [1x3]

for k_base = 1:num_baselines
    ant1_pos = all_baseline_ant_coords(k_base, 1:3);
    ant2_pos = all_baseline_ant_coords(k_base, 4:6);

    dist_to_ant1 = norm(S_guess_xyz - ant1_pos);
    dist_to_ant2 = norm(S_guess_xyz - ant2_pos);

    % 理论 DTOA (ns)
    theoretical_dtoa_ns = (dist_to_ant1 - dist_to_ant2) / c_mpns;

    errors(k_base) = (theoretical_dtoa_ns - all_measured_dtoas_ns(k_base)) / all_sigmas_ns(k_base);
end
end



function pulse_catalog = create_pulse_catalog(waveform, waveform2,waveform3,fs_hz, threshold, snippet_len)
% create_pulse_catalog: 从波形中检测所有脉冲并提取其特征。
%
% 输入:
%   waveform        - ch1通道的完整信号波形
%   waveform2        - ch2通道的完整信号波形
%   waveform3       - ch3通道的完整信号波形
%   fs_hz - 采率 (Hz), 例如 200e6
%   threshold       - findpeaks的振幅阈值
%   snippet_len     - 提取的波形片段长度 (点数), 例如 512
%
% 输出:
%   pulse_catalog   - 结构体数组, 每个元素代表一个脉冲及其特征


% --- Step 1: 使用 findpeaks 检测所有脉冲 ---
% 这里的MinPeakDistance需要根据您的密集程度仔细调整
[pks, locs] = findpeaks(waveform, 'MinPeakHeight', threshold, 'MinPeakDistance', snippet_len/4);

num_pulses = numel(locs);
if num_pulses == 0
    pulse_catalog = [];
    return;
end

% --- Step 2: 为每个脉冲提取特征 ---

% 初始化结构体数组
pulse_catalog = repmat(struct('loc', 0, 'amp', 0, 'snippet', [], 'snippet2', [],'snippet3', [],'fwhm_ns', 0, 'risetime_ns', 0), num_pulses, 1);

ts_ns = 1 / fs_hz * 1e9; % 每个采样点的时间 (ns)

for k = 1:num_pulses
    pk_loc = locs(k);
    pk_amp = pks(k);

    pulse_catalog(k).loc = pk_loc;
    pulse_catalog(k).amp = pk_amp;

    % 2.1 提取波形片段 (Snippet)
    s_start = max(1, pk_loc - floor(snippet_len/2));
    s_end   = min(length(waveform), pk_loc + floor(snippet_len/2) - 1);
    pulse_catalog(k).snippet = waveform(s_start:s_end);
    pulse_catalog(k).snippet2 = waveform2(s_start:s_end);
    pulse_catalog(k).snippet3 = waveform3(s_start:s_end);
    % 2.2 提取半峰全宽 (FWHM - Full Width at Half Maximum)
    half_max = pk_amp / 2;
    % 向左寻找半高点
    idx_left = find(waveform(1:pk_loc) < half_max, 1, 'last');
    % 向右寻找半高点
    idx_right = find(waveform(pk_loc:end) < half_max, 1, 'first') + pk_loc - 1;
    if ~isempty(idx_left) && ~isempty(idx_right)
        pulse_catalog(k).fwhm_ns = (idx_right - idx_left) * ts_ns;
    end

    % 2.3 提取上升时间 (10% - 90%)
    amp_10 = 0.1 * pk_amp;
    amp_90 = 0.9 * pk_amp;
    idx_10 = find(waveform(1:pk_loc) < amp_10, 1, 'last');
    idx_90 = find(waveform(1:pk_loc) < amp_90, 1, 'last');
    if ~isempty(idx_10) && ~isempty(idx_90)
        pulse_catalog(k).risetime_ns = (idx_90 - idx_10) * ts_ns;
    end
end

end

function [best_match_pulse, min_cost] = find_best_match(yld_pulse, chj_catalog, weights)
% find_best_match: 在CHJ目录中为给定的YLD脉冲寻找最佳匹配.
%
% 输入:
%   yld_pulse           - 单个YLD脉冲的结构体
%   chj_catalog         - 完整的CHJ脉冲目录
%   weights             - 包含w1,w2,w3的代价函数权重结构体
%
% 输出:
%   best_match_pulse    - CHJ目录中匹配上的脉冲结构体 (找不到则为空)
%   min_cost            - 匹配的最小代价值

best_match_pulse = [];
min_cost = inf;

% --- Step 3.1: 特征相似度匹配 ---
for i = 1:numel(chj_catalog)
    chj_candidate = chj_catalog(i);

    % 计算各项代价
    % Cost 1: 波形相似度 (1 - 归一化互相关系数)
    % 确保两个snippet等长
    len1 = length(yld_pulse.snippet);
    len2 = length(chj_candidate.snippet);
    if len1 ~= len2
        min_len = min(len1, len2);
        r_corr_max = max(xcorr(yld_pulse.snippet(1:min_len), chj_candidate.snippet(1:min_len), 'normalized'));
    else
        r_corr_max = max(xcorr(yld_pulse.snippet, chj_candidate.snippet, 'normalized'));
    end
    %     损失函数，损失，越小越好，而互相关是越大越好，可以通过1- 来转化为越小越好
    cost1 = 1 - r_corr_max;

    % Cost 2: 半峰全宽差异
    cost2 = abs(yld_pulse.fwhm_ns - chj_candidate.fwhm_ns);

    % Cost 3: 上升时间差异
    cost3 = abs(yld_pulse.risetime_ns - chj_candidate.risetime_ns);

    % 计算总代价
    total_cost = weights.w1 * cost1 + weights.w2 * cost2 + weights.w3 * cost3;

    if total_cost < min_cost
        min_cost = total_cost;
        best_match_pulse = chj_candidate;
    end
end
end


function [t12,t13,t23, Az_deg, El_deg, Rcorr, t123]  = get_2d_from_pulse(current_pulse,type)

Az_deg = 0;
El_deg = 0;
Rcorr = 0;
t123 = 100;
c = 0.299792458;
ch1 = current_pulse.snippet;
ch2 = current_pulse.snippet2;
ch3 = current_pulse.snippet3;
[ch1_up, ch2_up, ch3_up] = deal(...
    upsampling(ch1, 50)', ...
    upsampling(ch2, 50)', ...
    upsampling(ch3, 50)');
ch1_upsp = ch1_up(:,2);
ch2_upsp = ch2_up(:,2);
ch3_upsp = ch3_up(:,2);
%互相关
[r12_gcc,lags12_gcc] = xcorr(ch1_upsp,ch2_upsp,'normalized');
[r13_gcc,lags13_gcc] = xcorr(ch1_upsp,ch3_upsp,'normalized');
[r23_gcc,lags23_gcc] = xcorr(ch2_upsp,ch3_upsp,'normalized');
R12_gcc = max(r12_gcc);
R13_gcc = max(r13_gcc);
R23_gcc = max(r23_gcc);
t12_gcc = cal_tau(r12_gcc,lags12_gcc');
t13_gcc = cal_tau(r13_gcc,lags13_gcc');
t23_gcc = cal_tau(r23_gcc,lags23_gcc');
if strcmp(type, 'chj')
    % 从化局
    d12 = 41.6496;
    d13 = 36.9015;
    angle12 = -2.8381;
    angle13 = 50.3964;
    t12 = t12_gcc *0.1;
    t13 = t13_gcc *0.1+1.600061;
    t23 = t23_gcc *0.1+1.600061;
elseif strcmp(type, 'yld')
    % 引雷点
    angle12 = -110.8477;
    angle13 = -65.2405;
    d12 = 24.9586;
    d13 = 34.9335;
    t12 = t12_gcc *0.1;
    t13 = t13_gcc *0.1;
    t23 = t23_gcc *0.1;

end
cos_beta_0 =((c*t13*d12*sind(angle12))-(c*t12*sind(angle13)*d13))/(d13*d12*sind(angle12-angle13)) ;
cos_alpha_0 = ((c*t12)/d12-cos_beta_0*cosd(angle12))/sind(angle12);
if abs(cos_beta_0)>1 || abs(cos_alpha_0)>1
    return;
end
x0 = [cos_alpha_0,cos_beta_0];
% 调用lsqnonlin函数进行优化
options = optimoptions('lsqnonlin', 'MaxIter', 1000, 'TolFun', 1e-6);
x = lsqnonlin(@(x) objective(x, t12, t13, t23,type), x0, [-1 -1],[1 1], options);
% 输出最优的cos(α)和cos(β)值
cos_alpha_opt = x(1);
cos_beta_opt = x(2);
if abs(cos_alpha_opt)>1 || abs(cos_beta_opt)>1
    return;
end
Az = atan2( cos_alpha_opt,cos_beta_opt);
if abs(cos_beta_opt/cos(Az)) > 1
    return;
end
El = acos( cos_beta_opt/cos(Az) );
% 将弧度转换为角度
Az_deg = rad2deg(Az);
El_deg = rad2deg(El);
if Az_deg < 0
    Az_deg = Az_deg + 360;
end

t123 = t12 + t23 - t13;
Rcorr = (R12_gcc + R13_gcc + R23_gcc) / 3;
end

% 定义二维目标函数
function F = objective(x, t12_meas, t13_meas, t23_meas, type)
% 提取待优化的变量
cos_alpha = x(1);
cos_beta = x(2);

% 计算τij的理论值 τ_model (我将 obs 改为 model，语义更清晰)
tau_model = calculate_tau_obs(cos_alpha, cos_beta, type);

% t12, t13, t23 是测量的时延 (measurement)
% tau_model(1), tau_model(2), tau_model(3) 是根据当前 x 计算出的理论时延

% 计算残差向量
residual12 = t12_meas - tau_model(1);
residual13 = t13_meas - tau_model(2);
residual23 = t23_meas - tau_model(3);

% 返回残差向量 F
% lsqnonlin 会自动最小化 sum(F.^2)
F = [residual12; residual13; residual23];
end


% 定义计算τij的理想值τ_ij^obs的函数
function tau_ij_obs = calculate_tau_obs(cos_alpha, cos_beta, type)
% 初始化输出变量
tau_ij_obs = zeros(1, 3);

% 根据 type 参数选择不同的参数集
if strcmp(type, 'chj') % 从化局
    angle12 = -2.8381;
    angle13 = 50.3964;
    angle23 = 120.6568;
    d12 = 41.6496;
    d13 = 36.9015;
    d23 = 35.4481;
elseif strcmp(type, 'yld') % 引雷场
    angle12 = -110.8477;
    angle13 = -65.2405;
    angle23 = -19.6541;
    d12 = 24.9586;
    d13 = 34.9335;
    d23 = 24.9675;
else
    error('未知的类型：%s', type);
end

% 使用式(3)计算τij的理想值τ_ij^obs
tau_ij_obs(1) = (cos_alpha * sind(angle12) + cos_beta * cosd(angle12)) * d12 / 0.299792458;
tau_ij_obs(2) = (cos_alpha * sind(angle13) + cos_beta * cosd(angle13)) * d13 / 0.299792458;
tau_ij_obs(3) = (cos_alpha * sind(angle23) + cos_beta * cosd(angle23)) * d23 / 0.299792458;
end


function signal = read_signal(signal_path, r_length,r_loction)
fid  = fopen(signal_path,'r');%读取数据的位置

%使用fseek函数将文件指针移动到指定位置，以便读取数据。
%这里指定移动位置为r_location，表示移动到指定位置开始读取数据。
fseek(fid,r_loction*2,'bof');
%使用fread函数从文件中读取数据，读取的数据长度为r_length，数据以int16格式读取。
%将读取到的数据分别保存到变量ch_1、ch_2和ch_3中。
signal = fread(fid,r_length,'int16');
%关闭所有文件
fclose(fid);
end


%% 设计巴特沃斯带通滤波器
function filtered_signal = filter_bp(signal,f1,f2,order)
    Fs = 200e6;
    fn = Fs/2;
    Wn = [f1 f2]/fn;
    [b,a] = butter(order,Wn); 
    filtered_signal = filtfilt(b,a,signal);

end



%函数：对主窗口进行上采样
function new_signal = upsampling(original_signal,upsampling_factor)

% 原信号
original_x = (1:numel(original_signal))';
original_y = original_signal;
% 上采样后的采样点数
upsampled_length = length(original_x) * upsampling_factor;
% 上采样后的采样点的 x 坐标
upsampled_x = linspace(1, length(original_x), upsampled_length);
% 使用多项式插值对原信号进行上采样
interpolated_signal = interp1(original_x, original_y, upsampled_x, 'spline');
new_signal = [upsampled_x; interpolated_signal];
end



function tau = cal_tau(R, lag)
% 从数据中找到y的最大值及其索引
[~, max_index] = max(R);
tau = lag(max_index,1);
end
