%% ========================================================================
%  主程序: 脉冲级特征匹配定位算法
%  ========================================================================

clear; clc; close all;

%% --- 0. 初始化和参数定义 ---
% ... (引入您的所有基本变量: 站点位置, c, baseline_definitions等) ...
yld_sit = [0, 0, 0];
chj_sit = [2003.7972, -7844.7836, -27];
c = 0.299792458; % m/ns
all_baseline_definitions = ...; % 您的基线定义

% 信号和算法参数
sampling_rate = 200e6; % 假设采样率 200 MHz
ts_ns = 1 / sampling_rate * 1e9;
threshold = ...;       % findpeaks阈值
snippet_len = 512;     % 脉冲片段长度

% 代价函数权重 (需要仔细调试!)
weights.w1 = 0.7; % 波形相似度权重
weights.w2 = 0.2; % FWHM 权重
weights.w3 = 0.1; % 上升时间权重
cost_threshold = 0.5; % 只有总代价小于此阈值的匹配才被接受

% DTOA时间搜索窗
time_window_ns = 1000; % 在理论DTOA两侧搜索 +/- 1000 ns
time_window_samples = round(time_window_ns / ts_ns);

%% --- 1. 独立脉冲编目 (整个程序只执行一次) ---
% 读取并滤波您的YLD和CHJ的**完整**信号段
% filtered_yld_signal1 = ...;
% filtered_chj_signal1 = ...;

% 创建脉冲目录
yld_catalog = create_pulse_catalog(filtered_yld_signal1, sampling_rate, threshold, snippet_len);
chj_catalog = create_pulse_catalog(filtered_chj_signal1, sampling_rate, threshold, snippet_len);

% 检查目录是否为空
if isempty(yld_catalog) || isempty(chj_catalog)
    error('一个或两个站的脉冲目录为空，无法继续。');
end

%% --- 2. 遍历YLD脉冲目录进行匹配和定位 ---
num_yld_pulses = numel(yld_catalog);
all_S_results = [];
all_match_results = [];

% 初始化一个粗略的位置估计，用于计算第一个脉冲的理论DTOA
last_successful_S = yld_sit; % 可以是YLD站本身，或任何一个合理的先验位置

h = waitbar(0, '正在进行脉冲级匹配与定位...');

for i = 1:num_yld_pulses
    waitbar(i/num_yld_pulses, h);
    
    current_yld_pulse = yld_catalog(i);
    
    % --- 2.1 计算粗略的理论DTOA作为时间约束 ---
    t_chj_rough = norm(last_successful_S - chj_sit) / c;
    t_yld_rough = norm(last_successful_S - yld_sit) / c;
    rough_dtoa_ns = t_yld_rough - t_chj_rough;
    rough_dtoa_samples = round(rough_dtoa_ns / ts_ns);
    
    % --- 2.2 寻找最佳匹配 ---
    [matched_chj_pulse, match_cost] = find_best_match(current_yld_pulse, chj_catalog, ...
        rough_dtoa_samples, time_window_samples, weights);
    
    % --- 2.3 如果找到一个好的匹配，则进行精处理 ---
    if ~isempty(matched_chj_pulse) && match_cost < cost_threshold
        
        fprintf('YLD脉冲 #%d (loc:%d) 匹配成功 CHJ脉冲(loc:%d) with cost=%.3f\n', ...
            i, current_yld_pulse.loc, matched_chj_pulse.loc, match_cost);
            
        % a. 获取YLD和CHJ的2D测向结果
        %    这需要您改造 get_2d_result_single_window 函数，
        %    使其接收脉冲的snippet而不是整个数据块
        %    [yld_az, yld_el] = get_2d_from_pulse(current_yld_pulse, 'yld');
        %    [chj_az, chj_el] = get_2d_from_pulse(matched_chj_pulse, 'chj');
        
        % b. 进行三维交汇定位，得到S_initial
        %    ... (您的三维定位代码) ...
        
        % c. 进行DTOA优化和卡方检验
        %    all_measured_dtoas_ns(1:3) = ... YLD站内时差
        %    all_measured_dtoas_ns(4:6) = ... CHJ站内时差
        %    all_measured_dtoas_ns(7) = (current_yld_pulse.loc - matched_chj_pulse.loc) * ts_ns; % **直接测量的站间DTOA**
        
        %    ... (您完整的lsqnonlin和卡方检验代码) ...
        
        % d. 如果定位成功且质量好
        %    S_optimized = ...
        %    all_S_results = [all_S_results; S_optimized];
        %    all_match_results = [all_match_results; ...];
        %    
        %    % **用新的定位结果更新先验位置，为下一次匹配提供更准的约束**
        %    last_successful_S = S_optimized;
    end
end
close(h);

%% --- 3. 结果可视化 ---
% plot_3d(all_S_results, all_match_results);




function pulse_catalog = create_pulse_catalog(waveform, sampling_rate_hz, threshold, snippet_len)
% create_pulse_catalog: 从波形中检测所有脉冲并提取其特征。
%
% 输入:
%   waveform        - 单个通道的完整信号波形
%   sampling_rate_hz - 采率 (Hz), 例如 200e6
%   threshold       - findpeaks的振幅阈值
%   snippet_len     - 提取的波形片段长度 (点数), 例如 512
%
% 输出:
%   pulse_catalog   - 结构体数组, 每个元素代表一个脉冲及其特征

    fprintf('正在创建脉冲目录...\n');
    
    % --- Step 1: 使用 findpeaks 检测所有脉冲 ---
    % 这里的MinPeakDistance需要根据您的密集程度仔细调整
    [pks, locs] = findpeaks(waveform, 'MinPeakHeight', threshold, 'MinPeakDistance', 100);
    
    num_pulses = numel(locs);
    if num_pulses == 0
        pulse_catalog = [];
        fprintf('未在该波形中找到脉冲。\n');
        return;
    end
    
    fprintf('检测到 %d 个脉冲，正在提取特征...\n', num_pulses);
    
    % --- Step 2: 为每个脉冲提取特征 ---
    
    % 初始化结构体数组
    pulse_catalog = repmat(struct('loc', 0, 'amp', 0, 'snippet', [], 'fwhm_ns', 0, 'risetime_ns', 0), num_pulses, 1);
    
    ts_ns = 1 / sampling_rate_hz * 1e9; % 每个采样点的时间 (ns)
    
    for k = 1:num_pulses
        pk_loc = locs(k);
        pk_amp = pks(k);
        
        pulse_catalog(k).loc = pk_loc;
        pulse_catalog(k).amp = pk_amp;
        
        % 2.1 提取波形片段 (Snippet)
        s_start = max(1, pk_loc - floor(snippet_len/2));
        s_end   = min(length(waveform), pk_loc + floor(snippet_len/2) - 1);
        pulse_catalog(k).snippet = waveform(s_start:s_end);
        
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
    
    fprintf('脉冲目录创建完成。\n');
end


function [best_match_pulse, min_cost] = find_best_match(yld_pulse, chj_catalog, ...
    rough_dtoa_samples, time_window_samples, weights)
% find_best_match: 在CHJ目录中为给定的YLD脉冲寻找最佳匹配.
%
% 输入:
%   yld_pulse           - 单个YLD脉冲的结构体
%   chj_catalog         - 完整的CHJ脉冲目录
%   rough_dtoa_samples  - 粗略的理论DTOA (单位: 采样点数)
%   time_window_samples - 在粗略DTOA两侧的搜索窗口大小 (采样点数)
%   weights             - 包含w1,w2,w3的代价函数权重结构体
%
% 输出:
%   best_match_pulse    - CHJ目录中匹配上的脉冲结构体 (找不到则为空)
%   min_cost            - 匹配的最小代价值

    best_match_pulse = [];
    min_cost = inf;

    % --- Step 3.1: 时间窗约束 ---
    % 计算CHJ脉冲应该出现的理论位置范围
    expected_loc = yld_pulse.loc + rough_dtoa_samples;
    search_min_loc = expected_loc - time_window_samples;
    search_max_loc = expected_loc + time_window_samples;
    
    % 筛选出落在时间窗内的候选CHJ脉冲
    all_chj_locs = [chj_catalog.loc];
    candidate_indices = find(all_chj_locs >= search_min_loc & all_chj_locs <= search_max_loc);
    
    if isempty(candidate_indices)
        return; % 时间窗内没有候选者
    end
    
    % --- Step 3.2: 特征相似度匹配 ---
    for i = 1:numel(candidate_indices)
        chj_candidate = chj_catalog(candidate_indices(i));
        
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