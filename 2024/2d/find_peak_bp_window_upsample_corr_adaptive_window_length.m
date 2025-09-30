%% ========================================================================
%  二维定位算法 [四档自适应窗口最终版]
% =========================================================================
clear; clc; close all;

N = 3;
c = 0.299792458;
fs = 200e6;
step = 1e6;
upsampling_factor = 50;
start_signal_loc = 3.8e8;
end_signal_loc = 4.2e8;
signal_length = end_signal_loc - start_signal_loc;
all_start_signal_loc = start_signal_loc:step:end_signal_loc;
% % 从化局
% angle12 = -2.8381; angle13 = 50.3964; angle23 = 120.6568;
% d12 = 41.6496; d13 = 36.9015; d23 = 35.4481;
%引雷点
d12 = 24.9586; d13 = 34.9335; d23 = 24.9675;
angle12 = -110.8477; angle13 = -65.2405; angle23 = -19.6541;

%引雷点阈值
noise = read_signal('..\\20240822165932.6610CH1.dat',1e8,1e8);
filtered_noise = filter_bp(noise,30e6,80e6,5);
threshold = mean(filtered_noise)+5*std(filtered_noise);
% %从化局阈值
% noise = read_signal('..\\2024 822 85933.651462CH1.dat',1e8,1e8);
% filtered_noise = filter_bp(noise,30e6,80e6,5);
% threshold = mean(filtered_noise)+5*std(filtered_noise);
% 打开一个文本文件用于写入运行结果
fileID = fopen('result_yld_window_ADAPTIVE_1e6_factor4.txt', 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','Win_Len','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');
for j = 1:numel(all_start_signal_loc)-1

    current_block_start = all_start_signal_loc(j);
    current_block_end = all_start_signal_loc(j+1);

    fprintf('>>>>>> 正在处理信号块: %d -- %d \n', current_block_start, current_block_end);

    % --- 1. 读取当前处理块的完整信号 ---
    %     引雷点
    ch1 = read_signal('..\\20240822165932.6610CH1.dat', step, current_block_start);
    ch2 = read_signal('..\\20240822165932.6610CH2.dat', step, current_block_start);
    ch3 = read_signal('..\\20240822165932.6610CH3.dat', step, current_block_start);
    %     从化局
    % ch1 = read_signal('..\\2024 822 85933.651462CH1.dat',step,current_block_start);
    % ch2 = read_signal('..\\2024 822 85933.651462CH2.dat',step,current_block_start);
    % ch3 = read_signal('..\\2024 822 85933.651462CH3.dat',step,current_block_start+215/5);
    filtered_signal1 = filter_bp(ch1, 30e6, 80e6, 5);
    filtered_signal2 = filter_bp(ch2, 30e6, 80e6, 5);
    filtered_signal3 = filter_bp(ch3, 30e6, 80e6, 5);
    scout_pulse_catalog = find_pulses_advanced(filtered_signal1, 3.545, fs, 4);
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

        %
        %         %从化局
        %         t12 = t12_gcc *0.1;
        %         t13 = t13_gcc *0.1+1.600061;
        %         t23 = t23_gcc *0.1+1.600061;


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




function delta_t = delta_t(tij,tij_obs)
    delta_t = tij - tij_obs;
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



function filtered_signal = filter_fft(sig,f1,f2)
    y=fft(sig);%傅里叶变换得到一个复数
    fs = 200e6;
    n = length(y);
    %创建一个长度与输入信号 y 相同的零向量 yy；
    yy=zeros(1,length(y));
    % 使用 for 循环遍历信号 y 的每个采样点（m 表示当前的采样点索引，从0到 N-1）；
    for m=1:n-1
    %     判断当前采样点对应的频率是否在 8Hz 到 15Hz 范围内，如果在该范围内，则将对应的 yy 值置为0，表示该频率的信号被滤除；
        if m*(fs/n)<f1 || m*(fs/n)>f2 %将奈奎斯特之后的频率也滤除点掉
            yy(m+1)=0;
        else
    %         如果当前采样点对应的频率不在 8Hz 到 15Hz 范围内，则将 yy 的值保持为原始信号 y 的值。
            yy(m+1)=y(m+1);
        end
    end %将频率为8Hz-15Hz的信号的幅值置0
    filtered_signal=ifft(yy)';
    
end


% 设计巴特沃斯带通滤波器
function filtered_signal = filtersignal(signal,f1,f2,order,fs)
     % 滤波器通带下边界频率f1 滤波器通带上边界频率f2  滤波器阶数order
     % 创建滤波器对象
     filter = designfilt('bandpassiir', 'FilterOrder', order, 'HalfPowerFrequency1', f1, 'HalfPowerFrequency2', f2, 'SampleRate', fs);
     filtered_signal = filtfilt(filter,signal);
end





%函数：遍历整个信号，找到微尺度窗口下相关系数大于0.8的窗口
function correlated_windows = find_correlated_windows(signal1, signal2, signal3, window_size, threshold, Fs, N)
    % 窗口数量
    num_windows = N - window_size + 1;  
    % 存储相关系数大于阈值的窗口
    correlated_windows = [];
    for i = 1:num_windows
        % 提取当前窗口的数据
        window1 = signal1(i:i+window_size-1);
        window2 = signal2(i:i+window_size-1);
        window3 = signal3(i:i+window_size-1);
        [~,R12,~] = gccphat(window1,window2, Fs);
        [~,R13,~] = gccphat(window1,window3, Fs);
        [~,R23,~] = gccphat(window2,window3, Fs);
        max_R12 = maxvalue(R12);
        max_R13 = maxvalue(R13);
        max_R23 = maxvalue(R23);
        
        % 如果相关系数大于阈值，将窗口添加到结果列表中
        if max_R12 > threshold && max_R13 > threshold && max_R23 > threshold
            correlated_windows = [correlated_windows; i];
        end
    end
end

function filtered_signal = rfi_filter(ori_signal,sub_signal_length)
    filtered_signal = [];
    subsignal_starts = 1:sub_signal_length/2:length(ori_signal);
    for i = 1:length(subsignal_starts)
        if subsignal_starts(i) + sub_signal_length - 1 > length(ori_signal)
            continue
        end
        subsignal = ori_signal(subsignal_starts(i):subsignal_starts(i)+sub_signal_length-1);
        windowed_signal = window_plus(sub_signal_length,subsignal);
        sub_filtered_signal = datafilter(windowed_signal);
        filtered_signal = [filtered_signal; sub_filtered_signal(sub_signal_length*0.25+1:sub_signal_length*0.75)];
    end
end


%函数：寻找信号的最大峰值
function peaks = find_max_peaks(signal,threshold)
    % 找到信号中的峰值
    [pks,locs] = findpeaks(signal);
    % 根据阈值筛选峰值
    selectedLocs = locs(pks > threshold);
    % 获取过滤后的每个峰值的x值
    
    peaks =selectedLocs;
end


%函数：寻找信号的峰值
function peaks = find_peaks(signal,threshold)
    % 找到信号中的峰值
    [pks,locs] = findpeaks(signal);
    % 根据阈值筛选峰值
    peaks = locs(pks > threshold);
end


function matched_peaks_x = match_peaks(peaks1,peaks2,peaks3)
    matched_peaks_x = []; % 存储匹配峰值的x值矩阵
    for i = 1:numel(peaks1)
        curr_peak1 = peaks1(i);
        % 检查peaks2和peaks3中是否存在与peaks1对应的峰值且x值的差不大于4
        idx_peak2 = find(abs(peaks2 - curr_peak1) <= 10);  % 获取peaks2中匹配峰值的索引
        idx_peak3 = find(abs(peaks3 - curr_peak1) <= 10);  % 获取peaks3中匹配峰值的索引
        % 检查是否找到了匹配的峰值
        if ~isempty(idx_peak2) && ~isempty(idx_peak3)
            matched_peaks_x = [matched_peaks_x; [curr_peak1, peaks2(idx_peak2(1)), peaks3(idx_peak3(1))]];% 添加匹配峰值的x值矩阵
        end
    end
end


function max_index = maxindex(vector)
    % 提取实部部分
    
    max_value = max(vector);
    % 找到最大值对应的索引
    max_index = find(vector == max_value);
end


function mswed_signal = msw_signal(signal , peak_x ,length)
      % 找到峰值的 x 值在信号中的索引
    left_idx = max(peak_x - length+1, 1);  % 确定左边界的索引
    right_idx = min(peak_x + length, 10240);  % 确定右边界的索引
    mswed_signal = signal(left_idx:right_idx);  % 提取以中心 x 值为中心的左右40个采样点

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


function tau = showfitted(data)
    % 从数据中找到y的最大值及其索引
    [~, max_index] = max(data(:, 2));
    % 获取最大值周围的3个点的索引
    if max_index > length(data(:,2))-3 || max_index <3
        tau = 20/0.299552816 + 1;
    else
        fit_range = [ -3,-2,-1, 0, 1, 2, 3] + max_index;
        % 获取10个点的索引和对应的值
        fit_indices = data(fit_range, 1);
        fit_values = data(fit_range, 2);
        % 进行抛物线拟合
        coefficients = polyfit(fit_indices, fit_values, 2);
        % 根据拟合结果计算拟合曲线上的点
        fit_indices_curve = linspace(min(fit_indices), max(fit_indices), 1000);
        fit_values_curve = polyval(coefficients, fit_indices_curve);
        % 绘制原始数据和拟合曲线
%         figure;
%         plot(data(:, 1), data(:, 2))
%         plot(data(:, 1), data(:, 2), 'b', fit_indices_curve, fit_values_curve, 'r--');
%         legend('原始数据', '拟合曲线');
%         xlabel('y的索引');
%         ylabel('y的值');
        [~, max_index_fit] = max(fit_values_curve);
        tau = fit_indices_curve(1,max_index_fit);
    end
    
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


%对互相关函数进行上采样
function upsampling_gcc = upsampling_gc(r,lag,upsampling_factor)

    % 上采样后的采样点数
    upsampled_length = length(lag) * upsampling_factor;
    % 上采样后的采样点的 x 坐标
    upsampled_x = linspace(-numel(r)/2, numel(r)/2, upsampled_length);
    % 使用多项式插值对原信号进行上采样
    interpolated_signal = interp1(lag, r, upsampled_x, 'spline');
    upsampling_gcc = [upsampled_x; interpolated_signal]';

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
