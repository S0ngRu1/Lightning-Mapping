clear;
% === 全局参数 ===
N = 3;
c = 0.299792458; % 光速 (m/ns) - 用于主程序的时延计算
fs = 200e6;
fp_start = 20e6; % 通带起始 (Hz)
fp_end = 80e6;   % 通带结束 (Hz)
upsampling_factor = 50;
min_peak_distance = 512;
processing_window_len = 2048;
window_type = 'hann'; 

% === 几何与文件参数 ===
signal_length = 8e7; % 总处理长度
r_loction = 5e8;     % 起始绝对位置
d12 = 24.9586;
d13 = 34.9335;
d23 = 24.9675;
angle12 = -110.8477;
angle13 = -65.2405;
angle23 = -19.6541;

% === 分段处理参数 ===
block_size = 2e7;    % 每次处理 2000万点
overlap = 5000;      % 重叠长度
step_size = block_size - overlap; % 每次前进的步长

% === 结果文件初始化 ===
% 文件名增加标识
file_name = 'results\20230718175104_result_yld_5e8_5.7e8_window_2048_512_阈值4倍标准差_去零飘_'+string(fp_start/1e6)+'_'+string(fp_end/1e6)+'_'+ window_type +"_with_error" +'.txt';
fileID = fopen(file_name, 'w');

% 【修改点1】表头增加 Err_Az 和 Err_El
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','peak','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123', 'Err_Az', 'Err_El');

% === 阈值计算 ===
noise = read_signal('20230718175104.9180CH3.dat', 1e5, 1e8);
filtered_noise = filter_bp(noise, fp_start, fp_end, 5);
min_peak_height = mean(filtered_noise) + 3 * std(filtered_noise);
% 预先计算噪声功率用于SNR估计
noise_power_est = mean(filtered_noise.^2); 
clear noise filtered_noise; 

% === 准备窗函数 ===
if strcmp(window_type,'hann')
    win = hann(processing_window_len);
elseif strcmp(window_type,'exp_hann')
    win = exp_hanning(processing_window_len);
elseif strcmp(window_type,'blackman')
    win = blackman(processing_window_len);
elseif strcmp(window_type,'gausswin')
    win = gausswin(processing_window_len, 8);
else
    win = [];
end
if ~isempty(win), win = win(:); end

% === 进度条 ===
h = waitbar(0, '开始分段处理...');

% === 分段循环 ===
current_read_offset = 0; 
while current_read_offset < signal_length
    % 1. 计算当前读取长度
    current_len = min(block_size, signal_length - current_read_offset);
    if current_len < processing_window_len
        break;
    end
    waitbar(current_read_offset / signal_length, h, ...
        sprintf('正在处理分段: %.1f%% (Offset: %d)', (current_read_offset/signal_length)*100, current_read_offset));
    
    % 2. 读取数据
    current_abs_loc = r_loction + current_read_offset;
    ch1 = read_signal('20230718175104.9180CH1.dat', current_len, current_abs_loc);
    ch2 = read_signal('20230718175104.9180CH2.dat', current_len, current_abs_loc);
    ch3 = read_signal('20230718175104.9180CH3.dat', current_len, current_abs_loc);
    
    % 3. 滤波
    processed_ch1_yld = filter_bp(detrend(ch1), fp_start, fp_end, 5);
    processed_ch2_yld = filter_bp(detrend(ch2), fp_start, fp_end, 5);
    processed_ch3_yld = filter_bp(detrend(ch3), fp_start, fp_end, 5);
    clear ch1 ch2 ch3;
    
    % 4. 寻找峰值
    [pks, locs] = findpeaks(processed_ch1_yld, ...
        'MinPeakHeight', min_peak_height, ...
        'MinPeakDistance', min_peak_distance);
    num_peaks = length(locs);
    
    % 5. 遍历处理峰值
    for i = 1:num_peaks
        center_loc = locs(i); 
        
        % 边界检查与重叠去重
        start_idx = center_loc - floor(processing_window_len / 2);
        end_idx = center_loc + ceil(processing_window_len / 2) - 1;
        
        if start_idx < 1 || end_idx > current_len
            continue;
        end
        
        is_last_block = (current_read_offset + current_len >= signal_length);
        if ~is_last_block && (center_loc > (current_len - overlap))
            continue;
        end
        
        % 截取信号
        segment_ch1 = processed_ch1_yld(start_idx:end_idx);
        segment_ch2 = processed_ch2_yld(start_idx:end_idx);
        segment_ch3 = processed_ch3_yld(start_idx:end_idx);
        
        % 加窗
        if isempty(win)
            windowed_segment_ch1 = segment_ch1;
            windowed_segment_ch2 = segment_ch2;
            windowed_segment_ch3 = segment_ch3;
        else
            windowed_segment_ch1 = segment_ch1 .* win;
            windowed_segment_ch2 = segment_ch2 .* win;
            windowed_segment_ch3 = segment_ch3 .* win;
        end
        
        % 上采样
        [ch1_up, ch2_up, ch3_up] = deal(...
            upsampling(windowed_segment_ch1, upsampling_factor)', ...
            upsampling(windowed_segment_ch2, upsampling_factor)', ...
            upsampling(windowed_segment_ch3, upsampling_factor)');
        ch1_upsp = ch1_up(:,2);
        ch2_upsp = ch2_up(:,2);
        ch3_upsp = ch3_up(:,2);
        
        % 互相关
        [r12_gcc, lags12_gcc] = xcorr(ch1_upsp, ch2_upsp, 'normalized');
        [r13_gcc, lags13_gcc] = xcorr(ch1_upsp, ch3_upsp, 'normalized');
        [r23_gcc, lags23_gcc] = xcorr(ch2_upsp, ch3_upsp, 'normalized');
        R12_gcc = max(r12_gcc);
        R13_gcc = max(r13_gcc);
        R23_gcc = max(r23_gcc);
        t12_gcc = cal_tau(r12_gcc, lags12_gcc');
        t13_gcc = cal_tau(r13_gcc, lags13_gcc');
        t23_gcc = cal_tau(r23_gcc, lags23_gcc');
        
        % 还原时延 (单位：纳秒 ns)
        t12 = t12_gcc * 5 / upsampling_factor;
        t13 = t13_gcc * 5 / upsampling_factor;
        t23 = t23_gcc * 5 / upsampling_factor;
        
        % 几何解算
        cos_beta_0 = ((c*t13*d12*sind(angle12))-(c*t12*sind(angle13)*d13))/(d13*d12*sind(angle12-angle13));
        cos_alpha_0 = ((c*t12)/d12-cos_beta_0*cosd(angle12))/sind(angle12);
        
        if abs(cos_beta_0) > 1 || abs(cos_alpha_0) > 1
            continue;
        end
        
        x0 = [cos_alpha_0, cos_beta_0];
        options = optimoptions('lsqnonlin', 'MaxIter', 1000, 'TolFun', 1e-6, 'Display', 'off');
        x = lsqnonlin(@(x) objective(x, t12, t13, t23,'yld'), x0, [-1 -1], [1 1], options);
        cos_alpha_opt = x(1);
        cos_beta_opt = x(2);
        
        if abs(cos_alpha_opt) > 1 || abs(cos_beta_opt) > 1
            continue;
        end
        
        Az = atan2(cos_alpha_opt, cos_beta_opt);
        if abs(cos_beta_opt/cos(Az)) > 1
            continue;
        end
        El = acos(cos_beta_opt/cos(Az));
        Az_deg = rad2deg(Az);
        El_deg = rad2deg(El);
        if Az_deg < 0
            Az_deg = Az_deg + 360;
        end
        
        t123 = t12 + t23 - t13;
        Rcorr = (R12_gcc + R13_gcc + R23_gcc) / 3;
        global_start_loc = r_loction + current_read_offset + start_idx;

        % =========================================================================
        % 【修改点2】误差计算模块 (Corrected for Units)
        % =========================================================================
        
        % 1. 物理常数准备 (SI单位)
        c_SI = 2.99792458e8;  % 光速 m/s
        fs_Hz = fs;           % 采样率 Hz
        T_sec = processing_window_len / fs_Hz; % 窗口长度 (秒)
        f1_Hz = fp_start;     % Hz
        f2_Hz = fp_end;       % Hz
        
        % 2. 计算 SNR
        % 使用预先计算的 noise_power_est
        sig_power_ch1 = mean(windowed_segment_ch1.^2);
        sig_power_ch2 = mean(windowed_segment_ch2.^2);
        sig_power_ch3 = mean(windowed_segment_ch3.^2);
        avg_sig_power = (sig_power_ch1 + sig_power_ch2 + sig_power_ch3) / 3;
        
        SNR_linear = avg_sig_power / noise_power_est;
        if SNR_linear < 0.1, SNR_linear = 0.1; end
        
        % 3. 计算时延误差 Sigma_tau (Eq 3) -> 结果单位：秒
        term1 = 3 / (8 * pi^2);
        term2 = (1 + 2 * SNR_linear) / (SNR_linear^2);
        term3 = 1 / (T_sec * (f2_Hz^3 - f1_Hz^3));
        sigma_tau_sq = term1 * term2 * term3;
        sigma_tau_sec = sqrt(sigma_tau_sq); % 单位：秒
        
        % 4. 计算角度误差 (Eq 5)
        % 【关键转换】将纳秒时延转为秒，以匹配光速 c_SI
        t12_sec = t12 * 1e-9; 
        t13_sec = t13 * 1e-9;
        
        sum_tau_sq = t12_sec^2 + t13_sec^2;
        if sum_tau_sq < 1e-20, sum_tau_sq = 1e-20; end
        
        % 平均基线长度 (米)
        d_avg = (d12 + d13 + d23) / 3;
        
        % 方位角误差 (弧度)
        sigma_az_rad = (sigma_tau_sec * sqrt(t13_sec^2 + t12_sec^2)) / sum_tau_sq;
        
        max_possible_delay_sq = (d_avg / c_SI)^2;
        denom_sq = max_possible_delay_sq - sum_tau_sq;
        
        if denom_sq < 1e-20
            % 这种情况发生在仰角接近0度（地平线）时，误差趋向无穷大
            sigma_el_rad = deg2rad(20); % 截断
        else
            % 化简后的稳定公式
            sigma_el_rad = sigma_tau_sec / sqrt(denom_sq);
        end
        
        % 转为角度
        sigma_az_deg = rad2deg(sigma_az_rad);
        sigma_el_deg = rad2deg(sigma_el_rad);
        
        % 截断异常大值
        if sigma_az_deg > 20, sigma_az_deg = 20; end
        if sigma_el_deg > 20, sigma_el_deg = 20; end
        
        % =========================================================================

        % 【修改点3】写入文件，包含最后两列误差
        fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
            global_start_loc, processing_window_len, t12, t13, t23, cos_alpha_opt, cos_beta_opt, ...
            Az_deg, El_deg, Rcorr, t123, sigma_az_deg, sigma_el_deg);
            
    end 
    
    current_read_offset = current_read_offset + step_size;
end 
fclose(fileID);
close(h);
fprintf('处理完成。\n');

% === 辅助函数区域 ===
function new_signal = upsampling(original_signal,upsampling_factor)
    original_x = (1:numel(original_signal))';
    original_y = original_signal;
    upsampled_length = length(original_x) * upsampling_factor;
    upsampled_x = linspace(1, length(original_x), upsampled_length);
    interpolated_signal = interp1(original_x, original_y, upsampled_x, 'spline');
    new_signal = [upsampled_x; interpolated_signal];
end

function signal = read_signal(signal_path, r_length,r_loction)
    fid  = fopen(signal_path,'r');
    fseek(fid,r_loction*2,'bof');
    signal = fread(fid,r_length,'int16');
    fclose(fid);
end

function tau = cal_tau(R, lag)
    [~, max_index] = max(R);
    tau = lag(max_index,1);
end

function F = objective(x, t12_meas, t13_meas, t23_meas, type)
    cos_alpha = x(1);
    cos_beta = x(2);
    tau_model = calculate_tau_obs(cos_alpha, cos_beta, type);
    residual12 = t12_meas - tau_model(1);
    residual13 = t13_meas - tau_model(2);
    residual23 = t23_meas - tau_model(3);
    F = [residual12; residual13; residual23];
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
    % 注意：这里使用的是光速的 m/ns 单位，以匹配 t_meas 的纳秒单位
    c_ns = 0.299792458; 
    tau_ij_obs(1) = (cos_alpha * sind(angle12) + cos_beta * cosd(angle12)) * d12 / c_ns;
    tau_ij_obs(2) = (cos_alpha * sind(angle13) + cos_beta * cosd(angle13)) * d13 / c_ns;
    tau_ij_obs(3) = (cos_alpha * sind(angle23) + cos_beta * cosd(angle23)) * d23 / c_ns;
end

function w = exp_hanning(n, alpha)
    if nargin < 2, alpha = 3; end
    k = 0:n-1;
    hann_win = 0.5 - 0.5 * cos(2*pi*k/(n-1));
    exp_factor = exp(-alpha * (abs(k - (n-1)/2) ./ ((n-1)/2)).^2);
    w = hann_win .* exp_factor;
    w = w / max(w);
end

function filtered_signal = filter_bp(signal,f1,f2,order)
    Fs = 200e6;
    fn = Fs/2;
    Wn = [f1 f2]/fn;
    [b,a] = butter(order,Wn);
    filtered_signal = filtfilt(b,a,signal);
end