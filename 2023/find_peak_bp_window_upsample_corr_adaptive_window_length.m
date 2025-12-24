%% ========================================================================
%  二维定位算法 [四档自适应窗口最终版]
% =========================================================================
clear; clc; close all;

N = 3;
c = 0.299792458;
fs = 200e6;
step = 1e6;
upsampling_factor = 10;
start_signal_loc = 3e8;
end_signal_loc = 6e8;
signal_length = end_signal_loc - start_signal_loc;
all_start_signal_loc = start_signal_loc:step:end_signal_loc;
%引雷点
d12 = 24.9586; d13 = 34.9335; d23 = 24.9675;
angle12 = -110.8477; angle13 = -65.2405; angle23 = -19.6541;
%引雷点阈值
noise = read_signal('20230718175104.9180CH1.dat',1e5,1e8);
filtered_noise = filter_bp(noise,30e6,80e6,5);
threshold = mean(filtered_noise)+5*std(filtered_noise);
% 打开一个文本文件用于写入运行结果
fileID = fopen('results\20230718175104_result_yld_window_ADAPTIVE_1e6_factor4.txt', 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','Win_Len','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');
for j = 1:numel(all_start_signal_loc)-1

    current_block_start = all_start_signal_loc(j);
    current_block_end = all_start_signal_loc(j+1);

    fprintf('>>>>>> 正在处理信号块: %d -- %d \n', current_block_start, current_block_end);

    % --- 1. 读取当前处理块的完整信号 ---
    %     引雷点
    ch1 = read_signal('20230718175104.9180CH1.dat', step, current_block_start);
    ch2 = read_signal('20230718175104.9180CH2.dat', step, current_block_start);
    ch3 = read_signal('20230718175104.9180CH3.dat', step, current_block_start);
    filtered_signal1 = filter_bp(ch1, 30e6, 80e6, 5);
    filtered_signal2 = filter_bp(ch2, 30e6, 80e6, 5);
    filtered_signal3 = filter_bp(ch3, 30e6, 80e6, 5);
    scout_pulse_catalog = find_pulses_advanced(filtered_signal1, 3.545, fs, 4, 10);
    pulse_count_in_chunk = numel(scout_pulse_catalog);
    % 调用决策函数，获得当前信号块应该使用的窗口长度
    dynamic_window_len = get_adaptive_window_length_4tier(pulse_count_in_chunk, step);
    fprintf('      本块密度: %d, 决策窗口: %d\n', pulse_count_in_chunk, dynamic_window_len);

    [peaks_in_block, locs_in_block] = findpeaks(filtered_signal1, 'MinPeakHeight', threshold, 'MinPeakDistance', dynamic_window_len/4);
    if isempty(locs_in_block)
        fprintf('      在本块内未找到有效脉冲，跳过。\n');
        continue;
    end

    num_peaks_in_block = numel(locs_in_block);
    h = waitbar(0, sprintf('正在处理块内 %d 个峰值...', num_peaks_in_block));

    for pi = 1:num_peaks_in_block
        waitbar(pi / num_peaks_in_block, h);
        idx = locs_in_block(pi);
        % 确保峰值不超出信号范围
        % 使用 dynamic_window_len 截取窗口信号
        win_start_idx = max(1, idx - floor(dynamic_window_len / 2) + 1);
        win_end_idx = min(step, idx + floor(dynamic_window_len / 2));

        % 截取窗口信号
        signal1 = filtered_signal1(win_start_idx:win_end_idx);
        signal2 = filtered_signal2(win_start_idx:win_end_idx);
        signal3 = filtered_signal3(win_start_idx:win_end_idx);
        % 去直流分量并应用窗函数
        [ch1_new, ch2_new, ch3_new] = deal(...
            real(windowsignal(detrend(signal1))), ...
            real(windowsignal(detrend(signal2))), ...
            real(windowsignal(detrend(signal3))));

        % 上采样
        [ch1_up, ch2_up, ch3_up] = deal(...
            upsampling(ch1_new, upsampling_factor)', ...
            upsampling(ch2_new, upsampling_factor)', ...
            upsampling(ch3_new, upsampling_factor)');
        ch1_upsp = ch1_up(:,2);
        ch2_upsp = ch2_up(:,2);
        ch3_upsp = ch3_up(:,2);

        %互相关
        [r12_gcc,lags12_gcc] = xcorr(ch1_upsp ,ch2_upsp,'normalized');
        [r13_gcc,lags13_gcc] = xcorr(ch1_upsp ,ch3_upsp,'normalized');
        [r23_gcc,lags23_gcc] = xcorr(ch2_upsp ,ch3_upsp,'normalized');
        R12_gcc = max(r12_gcc);
        R13_gcc = max(r13_gcc);
        R23_gcc = max(r23_gcc);
        t12_gcc = cal_tau(r12_gcc,lags12_gcc');
        t13_gcc = cal_tau(r13_gcc,lags13_gcc');
        t23_gcc = cal_tau(r23_gcc,lags23_gcc');

        %引雷场
        t12 = t12_gcc *0.1;
        t13 = t13_gcc *0.1;
        t23 = t23_gcc *0.1;

        cos_beta_0 =((c*t13*d12*sind(angle12))-(c*t12*sind(angle13)*d13))/(d13*d12*sind(angle12-angle13)) ;
        cos_alpha_0 = ((c*t12)/d12-cos_beta_0*cosd(angle12))/sind(angle12);
        if abs(cos_beta_0)>1 || abs(cos_alpha_0)>1
            continue;
        end
        x0 = [cos_alpha_0,cos_beta_0];
        % 调用lsqnonlin函数进行优化
        options = optimoptions('lsqnonlin', 'MaxIter', 1000, 'TolFun', 1e-6);
        x = lsqnonlin(@(x) objective(x, t12, t13, t23,'yld'), x0, [-1 -1],[1 1], options);
        % 输出最优的cos(α)和cos(β)值
        cos_alpha_opt = x(1);
        cos_beta_opt = x(2);
        if abs(cos_alpha_opt)>1 || abs(cos_beta_opt)>1
            continue;
        end
        Az = atan2( cos_alpha_opt,cos_beta_opt);
        if abs(cos_beta_opt/cos(Az)) > 1
            continue;
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
        absolute_loc = current_block_start + idx;
        % 写入计算后的数据
        fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
            absolute_loc,dynamic_window_len, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
    end
    close(h);
end
% 关闭文件
fclose(fileID);



function tau = cal_tau(R, lag)
    % 从数据中找到y的最大值及其索引
    [~, max_index] = max(R);
    tau = lag(max_index,1);
end

% function tau = cal_tau(R, lag) % 新的、高精度的版本
%     [~, max_idx] = max(R);
% 
%     % 确保峰值不在数组的边缘，否则无法取到3个点
%     if max_idx == 1 || max_idx == length(R)
%         tau = lag(max_idx);
%         return;
%     end
% 
%     % 提取峰值点 (y2) 和它左右相邻的两个点 (y1, y3)
%     y1 = R(max_idx - 1);
%     y2 = R(max_idx);
%     y3 = R(max_idx + 1);
% 
%     % 抛物线顶点横坐标的偏移量公式： p = (y1 - y3) / (2 * (y1 - 2*y2 + y3))
%     % p 是相对于中心点 max_idx 的亚采样偏移量
%     % 注意：要处理分母为0或非常小的情况，避免计算错误
%     denominator = 2 * (y1 - 2*y2 + y3);
%     if abs(denominator) < 1e-9
%         p = 0; % 如果分母太小（例如，平顶），则不进行偏移
%     else
%         p = (y1 - y3) / denominator;
%     end
%     
%     % 计算最终的精确时延
%     % lag是等差数列，可以直接用 p 乘以步长
%     time_step = lag(2) - lag(1);
%     tau = lag(max_idx) + p * time_step;
% end

function fitted_peak_x = fitpeak(data,peak_index)
if peak_index+10 < 10240 && peak_index-10 > 0
    fit_range = (peak_index + (-10:10))';
elseif peak_index+6 < 10240 && peak_index-6 > 0
    fit_range = (peak_index + (-6:6))';
elseif peak_index+2 < 10240 && peak_index-2 > 0
    fit_range = (peak_index + (-2:2))';
else
    fitted_peak_x = peak_index;
    return;
end
fit_values = data(fit_range);
coefficients = polyfit(fit_range, fit_values, 2);
fit_indices_curve = linspace(min(fit_range), max(fit_range), 1000);
fit_values_curve = polyval(coefficients, fit_indices_curve);
% 绘制原始数据和拟合曲线
% figure;
% plot(1:length(data),data)
% plot(1:length(data),data, 'b', fit_indices_curve, fit_values_curve, 'r--');
% legend('原始数据', '拟合曲线');
% xlabel('y的索引');
% ylabel('y的值');
[~, max_index_fit] = max(fit_values_curve);
fitted_peak_x = fit_indices_curve(1,max_index_fit);
end


% 定义计算τij的理想值τ_ij^obs的函数
function tau_ij_obs = calculate_tau_obs(cos_alpha, cos_beta, type)
    % 初始化输出变量
    tau_ij_obs = zeros(1, 3);

    % 根据 type 参数选择不同的参数集
    if strcmp(type, 'chj') % 从化局
        angle12 = -2.8381;
        angle13 = 28.2006;
        angle23 = 87.3358;
        d12 = 41.6496;
        d13 = 48.5209;
        d23 = 25.0182;
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



%% 设计巴特沃斯带通滤波器
function filtered_signal = filter_bp(signal,f1,f2,order)
    Fs = 200e6;
    fn = Fs/2;
    Wn = [f1 f2]/fn;
    [b,a] = butter(order,Wn); 
    filtered_signal = filtfilt(b,a,signal);

end


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




function windowed_signal = windowsignal(signal)
%     r_length = length(signal);
%    % 使用汉明窗
%    window = hamming(r_length);
%    % 对滤波后的信号应用窗函数
%    windowed_signal = signal .* window; % 信号与窗函数相乘
% 
    X = fft(signal);      %变换到频域加窗
    r_length = length(X);
    window = hamming(r_length);
%     得到的是频域信号
    X_windowed = X .* window;

% %     % 进行逆傅里叶变换得到时域信号
      windowed_signal = ifft(X_windowed);

end




function delay = cal_delay(R_xy)
    r_xy = ifft(R_xy);

% 找到主峰值位置
[~, max_idx] = max(r_xy);

% 计算估计的时间延迟
delay = max_idx / 200e3;

end


function window_len = get_adaptive_window_length_4tier(pulse_count,single_length)
    thresholds = [180, 265, 483] * single_length/2e4;
    T1 = thresholds(1);
    T2 = thresholds(2);
    T3 = thresholds(3);
    
    if pulse_count <= T1
        window_len = 4096;
    elseif pulse_count > T1 && pulse_count <= T2
        window_len = 2048;
    elseif pulse_count > T2 && pulse_count <= T3
        window_len = 1024;
    else % pulse_count > T3
        window_len = 512;
    end
end



function pulse_catalog = find_pulses_advanced(waveform, noise_std, sampling_rate_hz, detection_threshold_factor,merge_gap_samples)
%
% 输入:
%   waveform        - 需要进行脉冲发现的信号波形
%   noise_std        - 噪声水平
%   sampling_rate_hz - 信号的采样率 (Hz)
%   detection_threshold_factor - 包络检测阈值因
%
% 输出:
%   pulse_catalog   - 结构体数组, 每个元素代表一个脉冲，包含:
%                     .start_idx      - 脉冲的起始采样点索引
%                     .end_idx        - 脉冲的结束采样点索引
%                     .peak_loc       - 脉冲内部最高峰的索引
%                     .precise_time_ns - 抛物线拟合得到的亚采样点精度时刻 (ns)
    
     %% --- 1. 计算希尔伯特包络和阈值 ---
    envelope = abs(hilbert(waveform));
    mean_envelope = mean(envelope);
    detection_threshold = noise_std * detection_threshold_factor;
    
    %% --- 2. 初步检测脉冲区域 ---
    pulse_regions = bwlabel(envelope > detection_threshold);
    num_initial_regions = max(pulse_regions);
    
    if num_initial_regions < 2
        regions_to_process = pulse_regions;
    else
        % ************* 融合邻近的脉冲区域 *************
        fprintf('初步检测到 %d 个区域，开始进行邻近区域融合...\n', num_initial_regions);
        
        merged_regions = pulse_regions; 
        
        while true % 持续循环，直到某一轮没有任何融合发生
            
            was_merged_in_this_pass = false; % 本轮融合的标志
            
            unique_labels = unique(merged_regions);
            unique_labels = unique_labels(unique_labels > 0); % 去掉背景0
            
            if numel(unique_labels) < 2
                break; % 如果只剩一个或没有区域，则结束
            end

            % 获取每个区域的边界
            stats = regionprops(merged_regions', 'BoundingBox');
            region_boundaries = zeros(numel(stats), 2);
            for k = 1:numel(stats)
                bbox = stats(k).BoundingBox;
                % --- 修正1：边界计算的“差一”错误 ---
                region_boundaries(k, 1) = floor(bbox(1));
                region_boundaries(k, 2) = floor(bbox(1)) + floor(bbox(3)) - 1; % 正确的终点
            end
            
            % 检查相邻区域的间隙
            for k = 1:numel(unique_labels) - 1
                current_label = unique_labels(k);
                next_label = unique_labels(k+1);
                
                % 找到属于这两个标签的所有点的索引
                current_indices = find(merged_regions == current_label);
                next_indices = find(merged_regions == next_label);
                
                current_region_end = max(current_indices);
                next_region_start = min(next_indices);
                
                gap = next_region_start - current_region_end;
                
                if gap < merge_gap_samples
                    % 将下一个区域合并到当前区域
                    merged_regions(merged_regions == next_label) = current_label;
                    was_merged_in_this_pass = true;
                    break; % 发生了一次融合，立即跳出内层循环，重新开始外层while循环
                end
            end
            
            if ~was_merged_in_this_pass
                % 如果完整的一轮扫描都没有发生任何融合，则说明融合已完成
                break;
            end
        end % while true
        
        % 对融合后的区域重新编号，使其连续
        [~, ~, regions_to_process] = unique(merged_regions);
        regions_to_process = reshape(regions_to_process, size(pulse_regions));
        regions_to_process(pulse_regions==0) = 0; % 恢复背景
    end
    
    num_final_regions = max(regions_to_process);
    
    %% --- 3. 对最终的脉冲区域进行精确界定和正时 (不变) ---
    pulse_catalog = [];
    
    for k = 1:num_final_regions
        region_indices = find(regions_to_process == k);
        if isempty(region_indices), continue; end
        [~, max_idx_in_region] = max(envelope(region_indices));
        peak_loc = region_indices(max_idx_in_region);
        start_idx = find_pulse_boundary(envelope, region_indices(1), mean_envelope, 'backward');
        end_idx = find_pulse_boundary(envelope, region_indices(end), mean_envelope, 'forward');
        if isnan(start_idx) || isnan(end_idx), continue; end
        precise_time_ns = get_precise_timing(envelope, peak_loc, sampling_rate_hz);
        pulse.start_idx = start_idx;
        pulse.end_idx = end_idx;
        pulse.peak_loc = peak_loc;
        pulse.precise_time_ns = precise_time_ns;
        if isempty(pulse_catalog), pulse_catalog = pulse; else, pulse_catalog(end+1) = pulse; end
    end
end

% --- 辅助函数1: 寻找脉冲边界 ---
function boundary_idx = find_pulse_boundary(envelope, peak_loc, mean_envelope, direction)
    search_loc = peak_loc;
    boundary_idx = NaN;
    max_iter = length(envelope); % 防止无限循环
    iter_count = 0;

    while iter_count < max_iter
        if strcmpi(direction, 'backward')
            % 检查峰值之前的5个点
            check_indices = (search_loc - 5) : (search_loc - 1);
            if check_indices(1) < 1
                boundary_idx = 1; % 到达信号开头
                break;
            end
            
            if all(envelope(check_indices) < mean_envelope)
                boundary_idx = check_indices(1);
                break;
            else
                search_loc = search_loc - 1;
            end
        else % forward
            % 检查峰值之后的5个点
            check_indices = (search_loc + 1) : (search_loc + 5);
            if check_indices(end) > length(envelope)
                boundary_idx = length(envelope); % 到达信号结尾
                break;
            end
            
            if all(envelope(check_indices) < mean_envelope)
                boundary_idx = check_indices(end);
                break;
            else
                search_loc = search_loc + 1;
            end
        end
        iter_count = iter_count + 1;
    end
end

% --- 辅助函数2: 抛物线拟合正时 ---
function precise_time_ns = get_precise_timing(envelope, peak_loc, sampling_rate_hz)
    ts_ns = 1 / sampling_rate_hz * 1e9;

    % 检查边界，确保可以取到5个点
    if peak_loc <= 2 || peak_loc >= length(envelope) - 1
        % 如果点太靠边，无法取5个点，则直接返回峰值点的时刻
        precise_time_ns = (peak_loc - 1) * ts_ns;
        return;
    end
    
    % 取以最高峰为中心的5个包络点
    y_fit = envelope(peak_loc-2 : peak_loc+2);
    x_fit = (-2:2)';
    
    % 进行二次多项式(抛物线)拟合
    p_coeffs = polyfit(x_fit, y_fit, 2);
    
    % 抛物线顶点公式 x = -b / (2a)
    % p_coeffs = [a, b, c]
    sub_sample_offset = -p_coeffs(2) / (2 * p_coeffs(1));
    
    % 计算亚采样点精度的位置
    precise_loc = peak_loc + sub_sample_offset;
    
    % 转换为纳秒
    precise_time_ns = (precise_loc - 1) * ts_ns;
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
