%% ========================================================================
%  二维定位算法 ,没有上采样
%  优化内容：Parfor并行、时域加窗、批量写入、预计算坐标
%  新增优化：互相关峰值抛物线拟合 (亚采样精度)
% =========================================================================
clear; clc; close all;

%% === 1. 全局参数设置 ===
N = 3;
c = 0.299792458; % 光速 (m/ns)
fs = 200e6;
step = 1e4;      % 步长
upsampling_factor = 1;
% 信号范围
start_signal_loc = 3.92e8;
end_signal_loc = 3.92e8+9e6;

all_start_signal_loc = start_signal_loc:step:end_signal_loc;
% 引雷点几何参数
d12 = 24.9586; d13 = 34.9335; d23 = 24.9675;
angle12 = -110.8477; angle13 = -65.2405; angle23 = -19.6541;
% 结果文件名
output_filename = 'results\result_yld_ADAPTIVE_parabolic.txt';

%% === 2. 预计算与初始化 ===
% 预计算阈值 (只读一次噪声段)
noise = read_signal('..\\20240822165932.6610CH1.dat', 1e5, 1e8);
filtered_noise = filter_bp(noise, 30e6, 80e6, 5);
noise_power_est = mean(filtered_noise.^2); 
threshold = mean(filtered_noise) + 3 * std(filtered_noise);
clear noise filtered_noise;

% 打开文件写入表头
fileID = fopen(output_filename, 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','Win_len','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123', 'Err_Az', 'Err_El');

% 启动并行池
if isempty(gcp('nocreate'))
    parpool; 
end

% 预定义优化选项
optim_opts = optimoptions('lsqnonlin', 'MaxIter', 100, 'TolFun', 1e-4, ...
    'Display', 'off', 'Algorithm', 'trust-region-reflective');

fprintf('开始处理，共 %d 个数据块...\n', numel(all_start_signal_loc)-1);
total_timer = tic;

%% === 3. 主循环 (分块处理) ===
for j = 1:numel(all_start_signal_loc)-1
    current_block_start = all_start_signal_loc(j);
    
    % --- 读取数据 ---
    ch1 = read_signal('..\\20240822165932.6610CH1.dat', step, current_block_start);
    ch2 = read_signal('..\\20240822165932.6610CH2.dat', step, current_block_start);
    ch3 = read_signal('..\\20240822165932.6610CH3.dat', step, current_block_start);
    
    % --- 滤波 ---
    filtered_signal1 = filter_bp(ch1, 30e6, 80e6, 5);
    filtered_signal2 = filter_bp(ch2, 30e6, 80e6, 5);
    filtered_signal3 = filter_bp(ch3, 30e6, 80e6, 5);
    
    % --- 自适应窗口决策 ---
    % 注意：find_pulses_advanced 和 get_adaptive_window_length_4tier 需确保在路径中
    scout_pulse_catalog = find_pulses_advanced(filtered_signal1, 3.545, fs, 3, 5);
    pulse_count_in_chunk = numel(scout_pulse_catalog);
    dynamic_window_len = get_adaptive_window_length_4tier(pulse_count_in_chunk, step);
    
    % --- 寻找峰值 ---
    [~, locs_in_block] = findpeaks(filtered_signal1, 'MinPeakHeight', threshold, 'MinPeakDistance', dynamic_window_len/4);
    
    num_peaks = numel(locs_in_block);
    if num_peaks == 0
        continue;
    end
    
    fprintf('块 %d/%d: 密度 %d, 窗口 %d, 峰值数 %d\n', j, numel(all_start_signal_loc)-1, pulse_count_in_chunk, dynamic_window_len, num_peaks);
    
    % --- 预计算辅助变量 ---
    win_vec = hamming(dynamic_window_len); 
    x_orig = (1:dynamic_window_len)';
    x_up = linspace(1, dynamic_window_len, dynamic_window_len * upsampling_factor);
    
    block_results = nan(num_peaks, 13);
    
    % === 4. 并行循环处理峰值 ===
    parfor pi = 1:num_peaks
        idx = locs_in_block(pi);
        
        % 边界检查
        half_win = floor(dynamic_window_len / 2);
        win_start_idx = idx - half_win + 1;
        win_end_idx = idx + half_win;
        
        if win_start_idx < 1 || win_end_idx > step
            continue;
        end
        
        % 截取数据
        sig1 = filtered_signal1(win_start_idx:win_end_idx);
        sig2 = filtered_signal2(win_start_idx:win_end_idx);
        sig3 = filtered_signal3(win_start_idx:win_end_idx);
        
        % 去趋势 + 时域加窗 
        sig1 = detrend(sig1) .* win_vec;
        sig2 = detrend(sig2) .* win_vec;
        sig3 = detrend(sig3) .* win_vec;
        
        % 上采样 (Spline插值)
        s1_up = interp1(x_orig, sig1, x_up, 'spline');
        s2_up = interp1(x_orig, sig2, x_up, 'spline');
        s3_up = interp1(x_orig, sig3, x_up, 'spline');
        
        % 互相关
        [r12, lags] = xcorr(s1_up, s2_up, 'normalized');
        [r13, ~]    = xcorr(s1_up, s3_up, 'normalized');
        [r23, ~]    = xcorr(s2_up, s3_up, 'normalized');
        
        % 寻找整数索引峰值
        [R12_max, max_idx12] = max(r12);
        [R13_max, max_idx13] = max(r13);
        [R23_max, max_idx23] = max(r23);
        
        % =================================================================
        % 【新增优化】 抛物线拟合 (Parabolic Fit) 计算亚采样偏移量 delta
        % =================================================================
        delta12 = calc_parabolic_offset(r12, max_idx12);
        delta13 = calc_parabolic_offset(r13, max_idx13);
        delta23 = calc_parabolic_offset(r23, max_idx23);
        
        % 纳秒时延计算 (包含小数部分 delta)
        % time_step = 1 / (fs * upsampling_factor) * 1e9 = 5ns / 50 = 0.1ns
        time_step = 5 / upsampling_factor; 
        
        t12 = (lags(max_idx12) + delta12) * time_step;
        t13 = (lags(max_idx13) + delta13) * time_step;
        t23 = (lags(max_idx23) + delta23) * time_step;
        % =================================================================
        
        % 初值计算 
        cb0 = ((c*t13*d12*sind(angle12))-(c*t12*sind(angle13)*d13))/(d13*d12*sind(angle12-angle13));
        ca0 = ((c*t12)/d12 - cb0*cosd(angle12))/sind(angle12);
        
        if abs(cb0) > 1 || abs(ca0) > 1
            continue; 
        end
        
        % 非线性优化
        x0 = [ca0, cb0];
        try
            x_opt = lsqnonlin(@(x) objective_fast(x, t12, t13, t23), x0, [-1 -1], [1 1], optim_opts);
        catch
            continue;
        end
        
        ca_opt = x_opt(1);
        cb_opt = x_opt(2);
        
        % 角度转换
        Az_rad = atan2(ca_opt, cb_opt);
        if abs(cb_opt/cos(Az_rad)) > 1
            continue;
        end
        El_rad = acos(cb_opt/cos(Az_rad));
        
        Az = rad2deg(Az_rad);
        El = rad2deg(El_rad);
        if Az < 0, Az = Az + 360; end
        
        % 误差计算
        t123 = t12 + t23 - t13;
        Rcorr_avg = (R12_max + R13_max + R23_max) / 3;
        
        % 理论误差估计
        sig_pow = (mean(sig1.^2) + mean(sig2.^2) + mean(sig3.^2)) / 3;
        snr_lin = max(0.1, sig_pow / noise_power_est);
        T_sec = dynamic_window_len / fs;
        term1 = 3 / (8 * pi^2);
        term2 = (1 + 2 * snr_lin) / (snr_lin^2);
        term3 = 1 / (T_sec * (80e6^3 - 30e6^3));
        sigma_tau = sqrt(term1 * term2 * term3);
        
        t12_s = t12 * 1e-9; t13_s = t13 * 1e-9;
        sum_tau_sq = max(1e-20, t12_s^2 + t13_s^2);
        err_az = rad2deg( (sigma_tau * sqrt(sum_tau_sq)) / sum_tau_sq );
        
        d_avg = (d12+d13+d23)/3; c_si = 2.99792e8;
        denom_el = (d_avg/c_si)^2 - sum_tau_sq;
        if denom_el < 1e-20
            err_el = 20;
        else
            err_el = rad2deg( sigma_tau / sqrt(denom_el) );
        end
        
        block_results(pi, :) = [current_block_start + idx, dynamic_window_len, ...
            t12, t13, t23, ca_opt, cb_opt, Az, El, Rcorr_avg, t123, err_az, err_el];
    end
    
    % --- 批量写入文件 ---
    block_results(any(isnan(block_results), 2), :) = [];
    if ~isempty(block_results)
        fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', block_results');
    end
end
fclose(fileID);
toc(total_timer);
fprintf('全部处理完成。\n');

%% === 辅助函数 ===

% 【新增】抛物线拟合计算偏移量
function delta = calc_parabolic_offset(r, idx)
    % r: 互相关序列
    % idx: 最大值的索引
    % delta: 距离 idx 的小数偏移量 (-0.5 到 0.5)
    
    % 边界保护
    if idx <= 1 || idx >= length(r)
        delta = 0;
        return;
    end
    
    y1 = r(idx-1);
    y2 = r(idx);
    y3 = r(idx+1);
    
    % 抛物线顶点公式: delta = (y1 - y3) / (2 * (y1 - 2*y2 + y3))
    denom = 2 * (y1 - 2*y2 + y3);
    
    if abs(denom) < 1e-10 % 防止分母为0 (极其平坦的峰)
        delta = 0;
    else
        delta = (y1 - y3) / denom;
    end
end

function signal = read_signal(signal_path, r_length, r_location)
    fid = fopen(signal_path, 'r');
    if fid == -1, error('无法打开文件: %s', signal_path); end
    fseek(fid, r_location * 2, 'bof');
    signal = fread(fid, r_length, 'int16');
    fclose(fid);
end

function filtered_signal = filter_bp(signal, f1, f2, order)
    persistent b a
    if isempty(b)
        Fs = 200e6;
        fn = Fs/2;
        Wn = [f1 f2]/fn;
        [b, a] = butter(order, Wn);
    end
    filtered_signal = filtfilt(b, a, signal);
end

function F = objective_fast(x, t12m, t13m, t23m)
    angle12 = -110.8477; angle13 = -65.2405; angle23 = -19.6541;
    d12 = 24.9586; d13 = 34.9335; d23 = 24.9675;
    c_ns = 0.299792458;
    ca = x(1); cb = x(2);
    
    t1_calc = (ca * sind(angle12) + cb * cosd(angle12)) * d12 / c_ns;
    t2_calc = (ca * sind(angle13) + cb * cosd(angle13)) * d13 / c_ns;
    t3_calc = (ca * sind(angle23) + cb * cosd(angle23)) * d23 / c_ns;
    
    F = [t12m - t1_calc; t13m - t2_calc; t23m - t3_calc];
end