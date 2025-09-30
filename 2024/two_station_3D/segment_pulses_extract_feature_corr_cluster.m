%% =======================================================================
%               高级闪电定位代码 (特征聚类分组法)
% ========================================================================
%
% 流程:
% 1. 脉冲分割: 从连续信号中分割出每个独立的脉冲波形。
% 2. 特征提取: 为每个脉冲计算时域和频域特征。
% 3. 脉冲聚类: 使用DBSCAN算法将特征相似的脉冲分组，认为它们同源。
% 4. 分组定位: 对每个脉冲簇(group)进行模板匹配和TDOA定位。
%
%==========================================================================

clear;
clc;
close all;

%% 1. 初始化参数
N = 3;
c = 0.299792458; % 光速 (m/ns)
fs = 200e6;      % 原始采样率 (Hz)
upsampling_factor = 50; % 上采样倍数

% -------- 引雷场几何参数 --------
d12 = 24.9586; d13 = 34.9335; d23 = 24.9675;
angle12 = -110.8477; angle13 = -65.2405; angle23 = -19.6541;

% -------- 读取和滤波数据 --------
signal_length = 2e7;
r_loction = 3.7e8;
ch1 = read_signal('..\\20240822165932.6610CH1.dat', signal_length, r_loction);
ch2 = read_signal('..\\20240822165932.6610CH2.dat', signal_length, r_loction);
ch3 = read_signal('..\\20240822165932.6610CH3.dat', signal_length, r_loction);

filtered_signal1 = filter_bp(ch1, 20e6, 80e6, 5);
filtered_signal2 = filter_bp(ch2, 20e6, 80e6, 5);
filtered_signal3 = filter_bp(ch3, 20e6, 80e6, 5);

% -------- 打开结果文件 --------
fileID = fopen('result_yld_cluster_method.txt', 'w');
fprintf(fileID, '%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Cluster_ID', 'Pulse_Count', 't12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt', 'Azimuth', 'Elevation', 'Rcorr', 't123');

%% =======================================================================
%                       新方法核心流程
% ========================================================================

%% 2. 脉冲分割 (代替旧的 findpeaks)
% 目标: 从信号中分割出所有独立的脉冲波形段
fprintf('Step 1: 正在分割脉冲...\n');
pulse_window_length = 512; % 定义单个脉冲的窗口长度
pulses = segment_pulses(filtered_signal1, filtered_signal2, filtered_signal3, fs, pulse_window_length);
fprintf('完成，共分割出 %d 个脉冲。\n', length(pulses));

%% 3. 多维特征提取
% 目标: 为每个脉冲计算特征向量(身份证)
fprintf('Step 2: 正在提取多维特征...\n');
num_pulses = length(pulses);
feature_matrix = zeros(num_pulses, 4); % N个脉冲 x 4个特征
h = waitbar(0, '正在提取特征...');
for i = 1:num_pulses
    % 我们以CH1的波形为基准提取特征
    feature_matrix(i, :) = extract_features(pulses(i).ch1, fs);
    waitbar(i/num_pulses, h);
end
close(h);
fprintf('完成特征提取。\n');

%% 4. 脉冲聚类 (DBSCAN)
% 目标: 将特征相似的脉冲分为一组
fprintf('Step 3: 正在进行DBSCAN聚类...\n');
% a. 特征标准化 (非常重要!)
X_normalized = zscore(feature_matrix);

% b. DBSCAN聚类 (参数需要根据数据特性调试)
epsilon = 0.05;      % 邻域半径
min_samples = 2;   % 形成簇的最小样本数
labels = dbscan(X_normalized, epsilon, min_samples);

unique_labels = unique(labels);
num_clusters = sum(unique_labels > 0);
fprintf('完成聚类，共发现 %d 个物理簇。\n', num_clusters);

%% 5. 分组定位 (新的主循环)
% 目标: 对每个脉冲簇进行定位
fprintf('Step 4: 正在对每个簇进行分组定位...\n');
cluster_ids = unique_labels(unique_labels > 0); % 获取所有有效的簇ID (忽略噪声-1)

h = waitbar(0, '正在处理脉冲簇...');
for i = 1:length(cluster_ids)
    current_cluster_id = cluster_ids(i);
    
    % a. 获取当前簇的所有脉冲
    indices_in_cluster = find(labels == current_cluster_id);
    pulse_count = length(indices_in_cluster);

    % b. 为簇计算一个高精度的TDOA (核心步骤)
    % 该函数内部会进行模板匹配，并返回时延信息
    [t12, t13, t23, Rcorr] = calculate_cluster_tdoa(pulses, indices_in_cluster, upsampling_factor, fs);
    
    % c. 进行质量控制 (Quality Control)
    if Rcorr < 0.5 % 对模板的相关性进行筛选
        continue;
    end
    t123 = t12 + t23 - t13;
    if abs(t123) > 1.0 % 对闭合时延进行筛选
       continue;
    end
    
    % d. 最小二乘优化求解 (与您原代码相同)
    cos_beta_0 = ((c*t13*d12*sind(angle12))-(c*t12*sind(angle13)*d13))/(d13*d12*sind(angle12-angle13));
    cos_alpha_0 = ((c*t12)/d12 - cos_beta_0*cosd(angle12))/sind(angle12);
    
    if abs(cos_beta_0) > 1 || abs(cos_alpha_0) > 1
        continue;
    end
    
    x0 = [cos_alpha_0, cos_beta_0];
    options = optimoptions('lsqnonlin', 'Display', 'off', 'MaxIter', 1000, 'TolFun', 1e-6);
    % **注意**: objective函数也需要修改以返回向量，这里假设已修改
    x = lsqnonlin(@(x) objective_vector(x, t12, t13, t23, 'yld', d12, d13, d23, angle12, angle13, angle23), x0, [-1 -1], [1 1], options);
    
    cos_alpha_opt = x(1);
    cos_beta_opt = x(2);
    
    if abs(cos_alpha_opt) > 1 || abs(cos_beta_opt) > 1
        continue;
    end
    
    Az = atan2(cos_alpha_opt, cos_beta_opt);
    if abs(cos_beta_opt / cos(Az)) > 1
        continue;
    end
    El = acos(cos_beta_opt / cos(Az));
    
    Az_deg = rad2deg(Az);
    if Az_deg < 0, Az_deg = Az_deg + 360; end
    El_deg = rad2deg(El);
    
    % e. 写入结果
    fprintf(fileID, '%-15d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
        current_cluster_id, pulse_count, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr, t123);

    waitbar(i/length(cluster_ids), h);
end
close(h);
fclose(fileID);
fprintf('定位完成，结果已保存至 result_yld_cluster_method.txt\n');

%% =======================================================================
%                       辅助函数 (Helper Functions)
% ========================================================================

%% 脉冲分割函数
function pulses = segment_pulses(s1, s2, s3, fs, pulse_window)
    % 使用合成包络信号进行初步峰值检测
    env1 = abs(hilbert(s1));
    env2 = abs(hilbert(s2));
    env3 = abs(hilbert(s3));
    combined_signal = env1 + env2 + env3;
    
    % 使用动态阈值和非常小的峰距来找到所有可能的脉冲
    noise_level = median(combined_signal) / 0.6745;
    dynamic_threshold = 3 * noise_level; % 较低的阈值以捕获更多脉冲
    min_peak_dist = 50; % 远小于脉冲窗口，确保密集脉冲不被漏掉
    
    [~, locs] = findpeaks(combined_signal, 'MinPeakHeight', dynamic_threshold, 'MinPeakDistance', min_peak_dist);
    
    % 创建一个结构体数组来存储所有分割出的脉冲
    pulses = struct('ch1', [], 'ch2', [], 'ch3', [], 'loc', []);
    count = 0;
    half_win = pulse_window / 2;
    
    for i = 1:length(locs)
        idx = locs(i);
        % 确保窗口不越界
        if idx > half_win && idx <= length(s1) - half_win
            count = count + 1;
            pulses(count).ch1 = s1(idx - half_win : idx + half_win - 1);
            pulses(count).ch2 = s2(idx - half_win : idx + half_win - 1);
            pulses(count).ch3 = s3(idx - half_win : idx + half_win - 1);
            pulses(count).loc = idx / fs; % 记录脉冲的绝对时间
        end
    end
end

%% 特征提取函数
function features = extract_features(waveform, fs)
    % 确保输入是列向量
    waveform = waveform(:);
    
    % --- 1. 时域特征 ---
    [peak_val, peak_idx] = max(abs(waveform));
    
    % 上升时间 (10% to 90%)
    thresh10 = 0.1 * peak_val;
    thresh90 = 0.9 * peak_val;
    idx10 = find(abs(waveform(1:peak_idx)) >= thresh10, 1, 'first');
    idx90 = find(abs(waveform(1:peak_idx)) >= thresh90, 1, 'first');
    if isempty(idx10) || isempty(idx90), rise_time = 0; else, rise_time = (idx90 - idx10) / fs; end
    
    % 脉宽 (50% to 50%)
    thresh50 = 0.5 * peak_val;
    idx50_rise = find(abs(waveform(1:peak_idx)) >= thresh50, 1, 'first');
    idx50_fall = find(abs(waveform(peak_idx:end)) <= thresh50, 1, 'first') + peak_idx - 1;
    if isempty(idx50_rise) || isempty(idx50_fall), pulse_width = 0; else, pulse_width = (idx50_fall - idx50_rise) / fs; end

    % --- 2. 频域特征 ---
    win = hann(length(waveform));
    [Pxx, f] = pwelch(waveform, win, 0, [], fs);
    
    % 频谱矩心
    f_centroid = sum(f .* Pxx) / sum(Pxx);
    
    % 频谱带宽
    f_bw = sqrt(sum(((f - f_centroid).^2) .* Pxx) / sum(Pxx));
    
    features = [rise_time, pulse_width, f_centroid, f_bw];
end

%% 为簇计算TDOA的函数 (模板匹配法)
function [t12, t13, t23, Rcorr_avg] = calculate_cluster_tdoa(pulses, indices, upsampling_factor, ~)
    num_in_cluster = length(indices);
    pulse_len = length(pulses(1).ch1);
    
    % 初始化存储矩阵
    cluster_s1 = zeros(num_in_cluster, pulse_len);
    cluster_s2 = zeros(num_in_cluster, pulse_len);
    cluster_s3 = zeros(num_in_cluster, pulse_len);
    
    % 收集簇内所有脉冲
    for i = 1:num_in_cluster
        pulse_idx = indices(i);
        cluster_s1(i, :) = pulses(pulse_idx).ch1;
        cluster_s2(i, :) = pulses(pulse_idx).ch2;
        cluster_s3(i, :) = pulses(pulse_idx).ch3;
    end
    
    % 创建模板 (通过对齐和平均)
    template_s1 = mean(align_waveforms(cluster_s1), 1);
    template_s2 = mean(align_waveforms(cluster_s2), 1);
    template_s3 = mean(align_waveforms(cluster_s3), 1);
    
    % 上采样模板
    template_s1_up = resample(template_s1, upsampling_factor, 1);
    template_s2_up = resample(template_s2, upsampling_factor, 1);
    template_s3_up = resample(template_s3, upsampling_factor, 1);
    

    % 使用GCC-PHAT进行互相关 (更鲁棒)
    % 使用 (:) 将输入强制转换为列向量，以符合函数要求
    [r12, lags12] = xcorr(template_s1_up(:), template_s2_up(:), 'normalized');
    [r13, lags13] = xcorr(template_s1_up(:), template_s3_up(:), 'normalized');
    [r23, lags23] = xcorr(template_s2_up(:), template_s3_up(:), 'normalized');

    % 使用抛物线插值计算高精度时延
    t12_lags = cal_tau_parabolic(r12, lags12);
    t13_lags = cal_tau_parabolic(r13, lags13);
    t23_lags = cal_tau_parabolic(r23, lags23);
    
    t12 = t12_lags*0.1; % 直接得到ns
    t13 = t13_lags*0.1;
    t23 = t23_lags*0.1;
    
    Rcorr_avg = (max(r12) + max(r13) + max(r23)) / 3;
end

%% 波形对齐辅助函数
function aligned_waveforms = align_waveforms(waveforms)
    [num_waves, len_wave] = size(waveforms);
    aligned_waveforms = zeros(num_waves, len_wave);
    center_point = floor(len_wave / 2);
    
    for i = 1:num_waves
        [~, max_idx] = max(abs(waveforms(i, :)));
        shift = center_point - max_idx;
        aligned_waveforms(i, :) = circshift(waveforms(i, :), shift);
    end
end

%% 高精度时延计算 (抛物线插值)
function tau_sec = cal_tau_parabolic(r, lags)
    [~, max_idx] = max(r);
    if max_idx == 1 || max_idx == length(r)
        tau_sec = lags(max_idx); return;
    end
    y1 = r(max_idx - 1); y2 = r(max_idx); y3 = r(max_idx + 1);
    denominator = 2 * (y1 - 2*y2 + y3);
    if abs(denominator) < 1e-9, p = 0; else, p = (y1 - y3) / denominator; end
    time_step = lags(2) - lags(1);
    tau_sec = lags(max_idx) + p * time_step;
end

%% 优化目标函数 (返回向量)
function F = objective_vector(x, t12_meas, t13_meas, t23_meas, ~, d12, d13, d23, angle12, angle13, angle23)
    cos_alpha = x(1);
    cos_beta = x(2);
    c_mns = 0.299792458; % m/ns

    t12_model = (d12 / c_mns) * (cos_beta * cosd(angle12) + cos_alpha * sind(angle12));
    t13_model = (d13 / c_mns) * (cos_beta * cosd(angle13) + cos_alpha * sind(angle13));
    t23_model = (d23 / c_mns) * (cos_beta * cosd(angle23) + cos_alpha * sind(angle23));
    
    F = [t12_model - t12_meas;
         t13_model - t13_meas;
         t23_model - t23_meas];
end