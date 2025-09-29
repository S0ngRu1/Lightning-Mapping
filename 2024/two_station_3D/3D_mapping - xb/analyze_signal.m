% =========================================================================
%                雷电信号特征对比分析脚本
% =========================================================================
% 目的: 量化对比两段雷电信号的关键特征，以解释其“可定位性”的差异。
%
% 分析维度:
%   1. 脉冲密度: 脉冲数量、率、脉冲间隔时间(IPI)分布。
%   2. 脉冲强度: 脉冲峰值幅度分布、信噪比(SNR)。
%   3. 频域特性: 功率谱密度(PSD)分布，是否存在干扰。
%   4. 波形质量: 通道间的互相关系数(Rcorr)和闭合时延误差(t123)分布。
% =========================================================================

clear;
clc;
close all;

%% 1. --- 配置区域 ---
% 设置两段信号的参数

% --- 通用分析参数 ---
params.fs = 200e6;
params.filter_band_under = 30e6; 
params.filter_band_up = 80e6; 
params.filter_order = 5;
params.upsampling_factor = 10; % 分析时可适当降低，加快速度
params.min_peak_dist_samples = 50; % 脉冲检测的最小峰距
params.corr_analysis_pulse_count = 1000; % 为节省时间，只对部分脉冲进行相关性分析

% --- 信号段 1 ("定位不准"的信号) ---
segment1.name = '噪声段';
segment1.file_ch1 = '..\\20240822165932.6610CH1.dat'; 
segment1.file_ch2 = '..\\20240822165932.6610CH2.dat';
segment1.file_ch3 = '..\\20240822165932.6610CH3.dat';
segment1.r_location = 2e8; 
segment1.length = 4e7;       

% --- 信号段 2 ("定位精确"的信号) ---
segment2.name = '稀疏/定位精确段';
segment2.file_ch1 = '..\\20240822165932.6610CH1.dat'; 
segment2.file_ch2 = '..\\20240822165932.6610CH2.dat';
segment2.file_ch3 = '..\\20240822165932.6610CH3.dat';
segment2.r_location = 4.5e8; % 
segment2.length = 4e7;       

%% 2. --- 主处理流程 ---
segments = {segment1, segment2};
results = cell(1, 2);

for i = 1:length(segments)
    fprintf('正在分析: %s...\n', segments{i}.name);
    
    % a. 读取数据
    fprintf('  - 读取数据...\n');
    ch1 = read_signal(segments{i}.file_ch1, segments{i}.length, segments{i}.r_location);
    ch2 = read_signal(segments{i}.file_ch2, segments{i}.length, segments{i}.r_location);
    ch3 = read_signal(segments{i}.file_ch3, segments{i}.length, segments{i}.r_location);
    
    % b. 滤波 (建议使用零相位滤波)
%     fprintf('  - 信号滤波...\n');
%     filt_ch1 = filter_bp(ch1, params.filter_band_under ,params.filter_band_up,params.filter_order);
%     filt_ch2 = filter_bp(ch2, params.filter_band_under ,params.filter_band_up,params.filter_order);
%     filt_ch3 = filter_bp(ch3, params.filter_band_under ,params.filter_band_up,params.filter_order);
    
    % c. 调用核心分析函数
    fprintf('  - 提取特征...\n');
    results{i} = analyze_lightning_segment(ch1, ch2, ch3, params);
    fprintf('分析完成.\n\n');
end

%% 3. --- 结果展示 ---

% a. 表格对比
fprintf('==================================================================\n');
fprintf('                     信号特征量化对比结果\n');
fprintf('==================================================================\n');
fprintf('%-25s | %-20s | %-20s\n', '特征指标', segments{1}.name, segments{2}.name);
fprintf('------------------------------------------------------------------\n');
fprintf('%-25s | %-20.0f | %-20.0f\n', '检测到的脉冲总数', results{1}.pulse_count, results{2}.pulse_count);
fprintf('%-25s | %-20.2f | %-20.2f\n', '平均脉冲率 (个/ms)', results{1}.pulse_rate_kHz*1000, results{2}.pulse_rate_kHz*1000);
fprintf('%-25s | %-20.2f | %-20.2f\n', '中位脉冲间隔 (us)', results{1}.ipi_median_us, results{2}.ipi_median_us);
fprintf('%-25s | %-20.2f | %-20.2f\n', '平均信噪比 (dB)', results{1}.snr_db, results{2}.snr_db);
fprintf('%-25s | %-20.2f | %-20.2f\n', '中位峰值幅度 (arb.)', results{1}.peak_median, results{2}.peak_median);
fprintf('------------------------------------------------------------------\n');
fprintf('%-25s | %-20.3f | %-20.3f\n', '平均互相关系数 Rcorr', results{1}.Rcorr_mean, results{2}.Rcorr_mean);
fprintf('%-25s | %-20.3f | %-20.3f\n', '中位互相关系数 Rcorr', results{1}.Rcorr_median, results{2}.Rcorr_median);
fprintf('%-25s | %-20.3f | %-20.3f\n', '闭合误差|t123|均值(ns)', results{1}.t123_mean_abs, results{2}.t123_mean_abs);
fprintf('==================================================================\n');

% b. 绘图对比
figure('Name', '信号特征对比分析', 'Position', [100, 100, 1200, 800]);

% 图1: 脉冲间隔时间 (IPI) 直方图
subplot(2, 2, 1);
histogram(results{1}.ipi_all_us, 'BinWidth', 1, 'Normalization', 'probability', 'DisplayName', segments{1}.name);
hold on;
histogram(results{2}.ipi_all_us, 'BinWidth', 1, 'Normalization', 'probability', 'DisplayName', segments{2}.name);
title('脉冲间隔时间(IPI)分布对比');
xlabel('IPI (us)'); ylabel('概率');
legend;
set(gca, 'XScale', 'log'); % 对数坐标更能看清短间隔
grid on;

% 图2: 功率谱密度 (PSD) 对比
subplot(2, 2, 2);
plot(results{1}.psd_freq / 1e6, 10*log10(results{1}.psd_power), 'DisplayName', segments{1}.name);
hold on;
plot(results{2}.psd_freq / 1e6, 10*log10(results{2}.psd_power), 'DisplayName', segments{2}.name);
title('功率谱密度(PSD)对比 (CH1)');
xlabel('频率 (MHz)'); ylabel('功率/频率 (dB/Hz)');
legend;
grid on;

% 图3: 互相关系数 (Rcorr) 直方图
subplot(2, 2, 3);
histogram(results{1}.Rcorr_all, 'BinWidth', 0.02, 'Normalization', 'probability', 'DisplayName', segments{1}.name);
hold on;
histogram(results{2}.Rcorr_all, 'BinWidth', 0.02, 'Normalization', 'probability', 'DisplayName', segments{2}.name);
title('互相关系数(Rcorr)分布对比');
xlabel('Rcorr'); ylabel('概率');
legend;
xlim([0, 1]);
grid on;

% 图4: 闭合时延误差 (|t123|) 直方图
subplot(2, 2, 4);
histogram(abs(results{1}.t123_all_ns), 'BinWidth', 0.2, 'Normalization', 'probability', 'DisplayName', segments{1}.name);
hold on;
histogram(abs(results{2}.t123_all_ns), 'BinWidth', 0.2, 'Normalization', 'probability', 'DisplayName', segments{2}.name);
title('闭合时延误差|t123|分布对比');
xlabel('|t123| (ns)'); ylabel('概率');
legend;
xlim([0, 10]);
grid on;

%% 4. --- 核心分析函数 ---
function stats = analyze_lightning_segment(ch1, ch2, ch3, params)
    % 初始化输出结构体
    stats = struct();
    fs = params.fs;

    % --- 脉冲检测 ---
    env1 = abs(hilbert(ch1));
    env2 = abs(hilbert(ch2));
    env3 = abs(hilbert(ch3));
    combined_env = env1 + env2 + env3;
    
    noise_level = median(combined_env) / 0.6745;
    dynamic_threshold = 3 * noise_level;
    [peaks, locs] = findpeaks(combined_env, 'MinPeakHeight', dynamic_threshold, 'MinPeakDistance', params.min_peak_dist_samples);
    
    % --- 1. 脉冲密度与强度分析 ---
    stats.pulse_count = length(peaks);
    duration_s = length(ch1) / fs;
    stats.pulse_rate_kHz = stats.pulse_count / (duration_s * 1000);
    
    ipi_all = diff(locs) / fs;
    stats.ipi_all_us = ipi_all * 1e6;
    stats.ipi_mean_us = mean(stats.ipi_all_us);
    stats.ipi_median_us = median(stats.ipi_all_us);
    
    stats.peak_median = median(peaks);
    stats.snr_db = 20 * log10(mean(peaks) / noise_level);
    
    % --- 2. 频域分析 (以CH1为例) ---
    [pxx, f] = pwelch(ch1, hann(4096), 2048, 4096, fs);
    stats.psd_power = pxx;
    stats.psd_freq = f;
    
    % --- 3. 波形质量与可定位性分析 ---
    % 为节省时间，只抽样部分脉冲进行分析
    num_to_analyze = min(stats.pulse_count, params.corr_analysis_pulse_count);
    rand_indices = randperm(stats.pulse_count, num_to_analyze);
    
    Rcorr_all = zeros(num_to_analyze, 1);
    t123_all = zeros(num_to_analyze, 1);
    
    window_half = 512; % 用于互相关的小窗口
    
    for i = 1:num_to_analyze
        idx = locs(rand_indices(i));
        if idx <= window_half || idx > length(ch1) - window_half
            continue;
        end
        
        s1 = ch1(idx - window_half : idx + window_half - 1);
        s2 = ch2(idx - window_half : idx + window_half - 1);
        s3 = ch3(idx - window_half : idx + window_half - 1);
        
        % 这里使用简化的互相关流程
        fs_up = fs * params.upsampling_factor;
        s1_up = resample(s1, params.upsampling_factor, 1);
        s2_up = resample(s2, params.upsampling_factor, 1);
        s3_up = resample(s3, params.upsampling_factor, 1);
        
        [r12, lags12] = xcorr(s1_up, s2_up, 'normalized');
        [r13, lags13] = xcorr(s1_up, s3_up, 'normalized');
        [r23, lags23] = xcorr(s2_up, s3_up, 'normalized');

        t12 = (lags12(r12 == max(r12))) * (1/fs_up) * 1e9;
        t13 = (lags13(r13 == max(r13))) * (1/fs_up) * 1e9;
        t23 = (lags23(r23 == max(r23))) * (1/fs_up) * 1e9;

        Rcorr_all(i) = (max(r12) + max(r13) + max(r23)) / 3;
        t123_all(i) = t12(1) + t23(1) - t13(1);
    end
    
    stats.Rcorr_all = Rcorr_all(Rcorr_all > 0); % 移除未计算的
    stats.t123_all_ns = t123_all(t123_all ~= 0);
    
    stats.Rcorr_mean = mean(stats.Rcorr_all);
    stats.Rcorr_median = median(stats.Rcorr_all);
    stats.t123_mean_abs = mean(abs(stats.t123_all_ns));
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
    fclose('all');
end

