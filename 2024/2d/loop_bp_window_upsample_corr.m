clear;
clc;
close all;

%% === 1. 参数设置 ===
N = 3;
c = 0.299792458; % 光速 m/ns
fs = 200e6;      % 采样率 Hz
fp_start = 30e6; % 通带起始 Hz
fp_end = 80e6;   % 通带结束 Hz
upsampling_factor = 2; % 上采样倍数 (使用抛物线拟合后，上采样可以设为1或者较小的值如5)
window_type = 'hann';

% --- 引雷点几何参数 ---
d12 = 24.9586; d13 = 34.9335; d23 = 24.9675;
angle12 = -110.8477; angle13 = -65.2405; angle23 = -19.6541;

% --- 信号处理参数 ---
signal_length = 1800000;
r_loction = 469200000;
min_window_distance = 256;
processing_window_len = 1024;
windows = 1:min_window_distance:signal_length-processing_window_len+1;

% --- 输入输出文件 ---
file_name = 'results\20240822165932_loop_result_yld_4.692e8_1.8e6_window_1024_256_ErrCalc_upsample_50.txt';

%% === 2. 预计算与噪声估计 ===
% 假设前 1e5 个点包含噪声
noise_len = 1e5;
noise_sig = read_signal('..\\20240822165932.6610CH1.dat', noise_len, 1e5); % 读取一段前面的数据
filtered_noise = filter_bp(detrend(noise_sig), fp_start, fp_end, 5);
noise_power_est = var(filtered_noise); % 估计噪声方差 (功率)
if noise_power_est == 0
    noise_power_est = 1e-10; % 防止除以零
end
fprintf('噪声功率估计完成: %.4e\n', noise_power_est);
clear noise_sig filtered_noise;

% 读取主信号
fprintf('正在读取主信号...\n');
ch1 = read_signal('..\\20240822165932.6610CH1.dat', signal_length, r_loction);
ch2 = read_signal('..\\20240822165932.6610CH2.dat', signal_length, r_loction);
ch3 = read_signal('..\\20240822165932.6610CH3.dat', signal_length, r_loction);

% 预处理 (去趋势 + 滤波)
processed_ch1_yld = filter_bp(detrend(ch1), fp_start, fp_end, 5);
clear ch1
processed_ch2_yld = filter_bp(detrend(ch2), fp_start, fp_end, 5);
clear ch2
processed_ch3_yld = filter_bp(detrend(ch3), fp_start, fp_end, 5);
clear ch3

% 打开文件写入表头
fileID = fopen(file_name, 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','Win_len','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123', 'Err_Az', 'Err_El');

% 预计算窗函数
win = hann(processing_window_len);
win = win(:);

% 优化选项预定义
options = optimoptions('lsqnonlin', 'MaxIter', 100, 'TolFun', 1e-6, 'Display', 'off');

%% === 3. 主循环处理 ===
h = waitbar(0, '正在处理窗口...');
num_windows = numel(windows);

for wi = 1:num_windows
    if mod(wi, 100) == 0
        waitbar(wi / num_windows, h, sprintf('正在处理窗口 %d/%d', wi, num_windows));
    end
    
    start_idx = windows(wi);
    end_idx = windows(wi) + processing_window_len - 1;
    
    % 截取并加窗
    segment_ch1 = processed_ch1_yld(start_idx:end_idx);
    segment_ch2 = processed_ch2_yld(start_idx:end_idx);
    segment_ch3 = processed_ch3_yld(start_idx:end_idx);
    
    windowed_ch1 = segment_ch1 .* win;
    windowed_ch2 = segment_ch2 .* win;
    windowed_ch3 = segment_ch3 .* win;
    
    % 上采样 (如果 factor=1，则不改变长度，upsampling函数需适配)
    % 注意：使用了抛物线拟合后，上采样倍数可以设小一点，或者设为1
    [ch1_up, ch2_up, ch3_up] = deal(...
        upsampling(windowed_ch1, upsampling_factor)', ...
        upsampling(windowed_ch2, upsampling_factor)', ...
        upsampling(windowed_ch3, upsampling_factor)');
    
    ch1_upsp = ch1_up(:,2);
    ch2_upsp = ch2_up(:,2);
    ch3_upsp = ch3_up(:,2);
    
    % 互相关
    [r12_gcc, lags12] = xcorr(ch1_upsp, ch2_upsp, 'normalized');
    [r13_gcc, lags13] = xcorr(ch1_upsp, ch3_upsp, 'normalized');
    [r23_gcc, lags23] = xcorr(ch2_upsp, ch3_upsp, 'normalized');
    
    % 获取最大值和整数索引
    [R12_max, idx12] = max(r12_gcc);
    [R13_max, idx13] = max(r13_gcc);
    [R23_max, idx23] = max(r23_gcc);
    
    % === [关键修改] 抛物线拟合计算亚采样偏移 ===
    delta12 = calc_parabolic_offset(r12_gcc, idx12);
    delta13 = calc_parabolic_offset(r13_gcc, idx13);
    delta23 = calc_parabolic_offset(r23_gcc, idx23);
    
    % 计算精确时延 (Lag + Delta) * Ts
    % 原始采样周期为 5ns (200MHz)，上采样后周期为 5/upsampling_factor
    dt_ns = 5 / upsampling_factor;
    
    t12_ns = (lags12(idx12) + delta12) * dt_ns;
    t13_ns = (lags13(idx13) + delta13) * dt_ns;
    t23_ns = (lags23(idx23) + delta23) * dt_ns;
    
    % 引雷场时延定义
    t12 = t12_ns;
    t13 = t13_ns;
    t23 = t23_ns;
    
    % 2D 初值计算
    cos_beta_0 = ((c*t13*d12*sind(angle12)) - (c*t12*sind(angle13)*d13)) / (d13*d12*sind(angle12-angle13));
    cos_alpha_0 = ((c*t12)/d12 - cos_beta_0*cosd(angle12)) / sind(angle12);
    
    if abs(cos_beta_0) > 1 || abs(cos_alpha_0) > 1
        continue;
    end
    
    % 非线性优化
    x0 = [cos_alpha_0, cos_beta_0];
    try
        x_opt = lsqnonlin(@(x) objective(x, t12, t13, t23, 'yld'), x0, [-1 -1], [1 1], options);
    catch
        continue;
    end
    
    cos_alpha_opt = x_opt(1);
    cos_beta_opt = x_opt(2);
    
    if abs(cos_alpha_opt) > 1 || abs(cos_beta_opt) > 1
        continue;
    end
    
    Az_rad = atan2(cos_alpha_opt, cos_beta_opt);
    if abs(cos_beta_opt/cos(Az_rad)) > 1
        continue;
    end
    El_rad = acos(cos_beta_opt/cos(Az_rad));
    
    Az_deg = rad2deg(Az_rad);
    El_deg = rad2deg(El_rad);
    if Az_deg < 0
        Az_deg = Az_deg + 360;
    end
    
    % === [关键修改] 误差与不确定度计算 ===
    t123 = t12 + t23 - t13;
    Rcorr_avg = (R12_max + R13_max + R23_max) / 3;
    
    % 1. 计算当前窗口信号平均功率 (三个通道平均)
    % 注意：使用加窗后的信号计算能量更准确反映参与互相关的部分，或者去窗能量补偿
    sig_pow = (mean(windowed_ch1.^2) + mean(windowed_ch2.^2) + mean(windowed_ch3.^2)) / 3;
    
    % 2. 计算线性信噪比 (SNR)
    snr_lin = max(0.1, sig_pow / noise_power_est);
    
    % 3. CRLB 时延不确定度 (sigma_tau)
    % T_sec: 信号持续时间 (秒)
    T_sec = processing_window_len / fs; 
    
    term1 = 3 / (8 * pi^2);
    term2 = (1 + 2 * snr_lin) / (snr_lin^2);
    % 带宽项: B^3 - b^3 (Hz)
    term3 = 1 / (T_sec * (fp_end^3 - fp_start^3));
    
    sigma_tau = sqrt(term1 * term2 * term3); % 单位：秒
    
    % 4. 传递误差到方位角 (Err_Az)
    % t12, t13 需要转换为秒进行计算
    t12_s = t12 * 1e-9; 
    t13_s = t13 * 1e-9;
    
    sum_tau_sq = t12_s^2 + t13_s^2;
    if sum_tau_sq < 1e-20
        sum_tau_sq = 1e-20; % 防止除零
    end
    
    % 根据您提供的公式: err = deg( (sigma * sqrt(sum)) / sum ) = deg( sigma / sqrt(sum) )
    err_az = rad2deg( (sigma_tau * sqrt(sum_tau_sq)) / sum_tau_sq );
    
    % 5. 传递误差到仰角 (Err_El)
    d_avg = (d12 + d13 + d23) / 3;
    c_si = 2.99792458e8; % 光速 m/s
    
    denom_el = (d_avg / c_si)^2 - sum_tau_sq;
    if denom_el <= 0
        err_el = 999; % 物理上不合理的情况，给一个大误差标记
    else
        err_el = rad2deg( sigma_tau / sqrt(denom_el) );
    end
    
    % 写入文件 (增加 Err_Az, Err_El)
    fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
        r_loction+start_idx, processing_window_len, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr_avg, t123, err_az, err_el);
end

fclose(fileID);
close(h);
fprintf('处理完成！\n');

%% === 辅助函数定义 ===

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

function signal = read_signal(signal_path, r_length, r_loction)
    fid = fopen(signal_path, 'r');
    if fid == -1
        error(['无法打开文件: ', signal_path]);
    end
    status = fseek(fid, r_loction * 2, 'bof');
    if status == -1
        fclose(fid);
        error('fseek 失败，位置可能超出文件范围');
    end
    signal = fread(fid, r_length, 'int16');
    fclose(fid);
end

function filtered_signal = filter_bp(signal, f1, f2, order)
    % 使用 persistent 变量避免重复设计滤波器，提高效率
    persistent b a
    if isempty(b)
        Fs = 200e6;
        fn = Fs/2;
        Wn = [f1 f2]/fn;
        [b, a] = butter(order, Wn);
    end
    filtered_signal = filtfilt(b, a, signal);
end

function new_signal = upsampling(original_signal, upsampling_factor)
    if upsampling_factor == 1
        % 如果因子为1，直接返回适量格式的数据
        original_x = (1:numel(original_signal))';
        new_signal = [original_x, original_signal];
        return;
    end
    
    original_x = (1:numel(original_signal))';
    original_y = original_signal;
    upsampled_length = length(original_x) * upsampling_factor;
    upsampled_x = linspace(1, length(original_x), upsampled_length);
    interpolated_signal = interp1(original_x, original_y, upsampled_x, 'spline');
    new_signal = [upsampled_x; interpolated_signal];
end

function tau_ij_obs = calculate_tau_obs(cos_alpha, cos_beta, type)
    tau_ij_obs = zeros(1, 3);
    if strcmp(type, 'chj') 
        angle12 = -2.8381; angle13 = 50.3964; angle23 = 120.6568;
        d12 = 41.6496; d13 = 36.9015; d23 = 35.4481;
    elseif strcmp(type, 'yld') 
        angle12 = -110.8477; angle13 = -65.2405; angle23 = -19.6541;
        d12 = 24.9586; d13 = 34.9335; d23 = 24.9675;
    else
        error('未知的类型：%s', type);
    end
    c_speed = 0.299792458;
    tau_ij_obs(1) = (cos_alpha * sind(angle12) + cos_beta * cosd(angle12)) * d12 / c_speed;
    tau_ij_obs(2) = (cos_alpha * sind(angle13) + cos_beta * cosd(angle13)) * d13 / c_speed;
    tau_ij_obs(3) = (cos_alpha * sind(angle23) + cos_beta * cosd(angle23)) * d23 / c_speed;
end

function F = objective(x, t12_meas, t13_meas, t23_meas, type)
    cos_alpha = x(1);
    cos_beta = x(2);
    tau_model = calculate_tau_obs(cos_alpha, cos_beta, type);
    F = [t12_meas - tau_model(1); t13_meas - tau_model(2); t23_meas - tau_model(3)];
end