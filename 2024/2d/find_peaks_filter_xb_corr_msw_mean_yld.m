N = 3;
c = 0.299792458;
fs = 200e6;
upsampling_factor = 50;
window_length = 1024;
bigwindows_length = window_length+100;
window = window_length * upsampling_factor;
msw_length = 50;
%引雷点
signal_length = 8e5;
r_loction = 4.694e8;
d12 = 24.9586;
d13 = 34.9335;
d23 = 24.9675;
angle12 = -110.8477;
angle13 = -65.2405;
angle23 = -19.6541;
ch1 = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction);
ch2 = read_signal('..\\20240822165932.6610CH2.dat',signal_length,r_loction);
ch3 = read_signal('..\\20240822165932.6610CH3.dat',signal_length,r_loction);


filtered_signal1 = filter_xb(ch1);
filtered_signal2 = filter_xb(ch2);
filtered_signal3 = filter_xb(ch3);


noise = read_signal('..\\20240822165932.6610CH1.dat',60000,2e8);
filtered_noise = filter_xb(noise);
threshold = 0.5*std(filtered_noise);

% 打开一个文本文件用于写入运行结果
fileID = fopen('result_yld_4.5-5.5e8——xb.txt', 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','peak','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');


% % 设置动态阈值，每4000点取一个阈值
% subsignal_length = 4000;
% subsignal_start = 1:subsignal_length:length(filtered_signal1);
% for subi = 1:numel(subsignal_start)
% subsignal1 = filtered_signal1(subsignal_start(subi):subsignal_start(subi)+subsignal_length-1);
% % threshold = 0.5 * mean(abs(subsignal1));
% threshold =  mean(abs(subsignal1)) + 3*std(subsignal1);

% 寻找峰值
[peaks, locs] = findpeaks(filtered_signal1, 'MinPeakHeight', threshold, 'MinPeakDistance', window_length/4);

% 存储所有峰值和阈值
all_peaks = peaks;
all_locs = locs;

% 遍历所有峰值
num_peaks = numel(all_peaks);
% 创建进度条
h = waitbar(0, '正在处理峰值...');
for pi = 1:num_peaks
    waitbar(pi / num_peaks, h, sprintf('正在处理峰值 %d/%d', pi, num_peaks));
    idx = all_locs(pi);

    % 确保峰值不超出信号范围
    if idx - (bigwindows_length / 2 - 1) <= 0 || idx + (bigwindows_length / 2) > length(filtered_signal1)
        continue;
    end

    % 截取窗口信号
    [signal1, signal2, signal3] = deal(...
        filtered_signal1(idx-(bigwindows_length/2-1):idx+(bigwindows_length/2)), ...
        filtered_signal2(idx-(bigwindows_length/2-1):idx+(bigwindows_length/2)), ...
        filtered_signal3(idx-(bigwindows_length/2-1):idx+(bigwindows_length/2)));
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

    %取窗口做互相关
    ch1_new = ch1_upsp(bigwindows_length*upsampling_factor/2-(window/2-1):bigwindows_length*upsampling_factor/2+(window/2));
    ch2_new = ch2_upsp(bigwindows_length*upsampling_factor/2-(window/2-1):bigwindows_length*upsampling_factor/2+(window/2));
    ch3_new = ch3_upsp(bigwindows_length*upsampling_factor/2-(window/2-1):bigwindows_length*upsampling_factor/2+(window/2));
    %互相关
    [r12_gcc,lags12_gcc] = xcorr(ch1_new,ch2_new,'normalized');
    [r13_gcc,lags13_gcc] = xcorr(ch1_new,ch3_new,'normalized');
    [r23_gcc,lags23_gcc] = xcorr(ch2_new,ch3_new,'normalized');
    R12_gcc = max(r12_gcc);
    R13_gcc = max(r13_gcc);
    R23_gcc = max(r23_gcc);
    t12_gcc = cal_tau(r12_gcc,lags12_gcc');
    t13_gcc = cal_tau(r13_gcc,lags13_gcc');
    t23_gcc = cal_tau(r23_gcc,lags23_gcc');
    ismsw = 0;
    if R12_gcc > 0.8 && R13_gcc > 0.8
        %如果R大于某个值则将根据互相关曲线的最大值进行数据平移的结果
        shifted_ch1 = ch1_upsp(bigwindows_length*upsampling_factor/2-(window/2-1):bigwindows_length*upsampling_factor/2+(window/2));
        shifted_ch2 = shift_signal(ch2_upsp,t12_gcc);
        shifted_ch2 = shifted_ch2(bigwindows_length*upsampling_factor/2-(window/2-1):bigwindows_length*upsampling_factor/2+(window/2));
        shifted_ch3 = shift_signal(ch3_upsp,t13_gcc);
        shifted_ch3 = shifted_ch3(bigwindows_length*upsampling_factor/2-(window/2-1):bigwindows_length*upsampling_factor/2+(window/2));

        threshold1 = mean(shifted_ch1);
        threshold2 = mean(shifted_ch2);
        threshold3 = mean(shifted_ch3);
        peaks1= find_peaks(shifted_ch1,threshold1);
        peaks2= find_peaks(shifted_ch2,threshold2);
        peaks3= find_peaks(shifted_ch3,threshold3);

        matched_peaks_x = match_peaks(peaks1,peaks2,peaks3);
        R12s = []; %用于存储R12s值的空向量
        R13s = []; %用于存储R13s值的空向量
        R23s = []; %用于存储R23s值的空向量
        t12s = [];
        t13s = [];
        t23s = [];
        if size(matched_peaks_x,1) ~= 0
            for i = 1 : size(matched_peaks_x, 1)
                %微尺度
                ch1_msw = msw_signal(shifted_ch1 , matched_peaks_x(i,1) ,msw_length,window);
                ch2_msw = msw_signal(shifted_ch2 , matched_peaks_x(i,2) ,msw_length,window);
                ch3_msw = msw_signal(shifted_ch3 , matched_peaks_x(i,3) ,msw_length,window);
                if numel(ch1_msw)~= msw_length*2 || numel(ch2_msw)~= msw_length*2 || numel(ch3_msw)~= msw_length*2
                    continue;
                end
                %对微尺度进行互相关
                [R12_msw,lags12_msw] = xcorr(ch1_msw,ch2_msw,'normalized');
                [R13_msw,lags13_msw] = xcorr(ch1_msw,ch3_msw,'normalized');
                [R23_msw,lags23_msw] = xcorr(ch2_msw,ch3_msw,'normalized');
                if max(R12_msw) > 0.8 && max(R13_msw) > 0.8 && max(R23_msw) > 0.8
                    t12_msw = cal_tau(R12_msw,lags12_msw');
                    t13_msw = cal_tau(R13_msw,lags13_msw');
                    t23_msw = cal_tau(R23_msw,lags23_msw');
                    R12s = [R12s R12_msw];
                    R13s = [R13s R13_msw];
                    R23s = [R23s R23_msw];
                    t12s = [t12s t12_msw];
                    t13s = [t13s t13_msw];
                    t23s = [t23s t23_msw];
                end
            end
            if size(t12s,1)~=0 && size(t13s,1)~=0 && size(t23s,1)~=0
                %引雷场
                t12 = (t12_gcc + mean(t12s))*0.1;
                t13 = (t13_gcc + mean(t13s))*0.1;
                t23 = (t23_gcc + mean(t23s))*0.1;

                cos_beta_0 =((c*t13*d12*sind(angle12))-(c*t12*sind(angle13)*d13))/(d13*d12*sind(angle12-angle13)) ;
                cos_alpha_0 = ((c*t12)/d12-cos_beta_0*cosd(angle12))/sind(angle12);
                if abs(cos_beta_0)>1 || abs(cos_alpha_0)>1
                    continue;
                end

                x0 = [cos_alpha_0,cos_beta_0];
                % 调用lsqnonlin函数进行优化
                options = optimoptions('lsqnonlin', 'MaxIter', 1000, 'TolFun', 1e-6);
                x = lsqnonlin(@(x) objective(x, t12, t13, t23), x0, [-1 -1],[1 1], options);
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

                % 写入计算后的数据
                fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
                     r_loction+idx-window/100,window/100, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
                ismsw = ismsw + 1;
            end
        end
    end
    if ismsw == 0
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
        x = lsqnonlin(@(x) objective(x, t12, t13, t23), x0, [-1 -1],[1 1], options);
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

        % 写入计算后的数据
        fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
             r_loction+idx-window/100,window/100, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
    end
end
% 关闭文件
fclose(fileID);
% 关闭进度条
close(h);


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


% 定义目标函数
function F = objective(x,t12,t13,t23)
    % 提取待优化的变量
    cos_alpha = x(1);
    cos_beta = x(2);

    % 计算τij的理想值τ_ij^obs
    tau_ij_obs = calculate_tau_obs(cos_alpha, cos_beta);
    % 计算Δt12, Δt13, Δt23
    delta_t12 = delta_t(t12,tau_ij_obs(1));
    delta_t13 = delta_t(t13,tau_ij_obs(2));
    delta_t23 = delta_t(t23,tau_ij_obs(3));

    % 计算目标函数，即式(4)
    F = (delta_t12^2 + delta_t13^2 + delta_t23^2) / 75;
end

function tau_ij_obs = calculate_tau_obs(cos_alpha, cos_beta)
    % 初始化输出变量
    tau_ij_obs = zeros(1, 3);
        % 引雷场
        angle12 = -110.8477;
        angle13 = -65.2405;
        angle23 = -19.6541;
        d12 = 24.9586;
        d13 = 34.9335;
        d23 = 24.9675;

    % 使用式(3)计算τij的理想值τ_ij^obs
    tau_ij_obs(1) = (cos_alpha * sind(angle12) + cos_beta * cosd(angle12)) * d12 / 0.299792458;
    tau_ij_obs(2) = (cos_alpha * sind(angle13) + cos_beta * cosd(angle13)) * d13 / 0.299792458;
    tau_ij_obs(3) = (cos_alpha * sind(angle23) + cos_beta * cosd(angle23)) * d23 / 0.299792458;
end



function tau = cal_tau(R, lag)
    % 从数据中找到y的最大值及其索引
    [~, max_index] = max(R);
    tau = lag(max_index,1);
end

function delta_t = delta_t(tij,tij_obs)
    delta_t = tij - tij_obs;
end

function filtered_signal = filter_xb(signal)
% 通用小波
    level = 4;          % 分解层数
    wavelet = 'db4';    % 小波基
    [c4, l4] = wavedec(signal, level, wavelet);
    sigma = median(abs(c4))/0.6745;
%     thr = sigma*sqrt(2*log(length(signal)));  % 通用阈值
    thr = sigma*sqrt(log(length(signal)));
    c_denoised4 = wthresh(c4, 's', 1.5*thr);  % 使用小波软阈值
    filtered_signal = normalize(waverec(c_denoised4, l4, wavelet));


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



function mswed_signal = msw_signal(signal , peak_x ,length)
      % 找到峰值的 x 值在信号中的索引
    left_idx = max(peak_x - length+1, 1);  % 确定左边界的索引
    right_idx = min(peak_x + length, 10240);  % 确定右边界的索引
    mswed_signal = signal(left_idx:right_idx);  % 提取以中心 x 值为中心的左右40个采样点

end



function shifted_signal = shift_signal(signal, shift_amount)

    % 使用 circshift 进行平移
    shifted_signal = circshift(signal, shift_amount);
    % 如果是向左平移，右侧补零；如果是向右平移，左侧补零
    if shift_amount < 0
        shifted_signal(end+shift_amount+1:end) = 0;
    else
        shifted_signal(1:shift_amount) = 0;
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
