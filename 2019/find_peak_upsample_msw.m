%读取数据
signal_length = 1e8;
ch1 = read_signal('../cross-correlation/20190604164852.7960CH1.dat',signal_length);
ch2 = read_signal('../cross-correlation/20190604164852.7960CH2.dat',signal_length);
ch3 = read_signal('../cross-correlation/20190604164852.7960CH3.dat',signal_length);

% x = downsample(filtered_signal1,50);
% filtered_signalx = preprocess(x);

N = 3;
d = 20;
c = 0.299792458;
window_length = 1024;

angle12 = -150;
angle13 = -90;
angle23 = -30;

% 打开一个文本文件用于写入运行结果
fileID = fopen('result14.txt', 'w');
% 写入第一行的数据介绍
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
     'Start_loc','Peak1_loc','Peak2_loc','Peak3_loc','t12', 't13', 't23', 'cosα', 'cosβ', 'Azimuth', 'Elevation', 'Rcorr', 't123');


        filtered_signal1 = filter_bp(ch1, 20e6 ,80e6 ,8);
        filtered_signal2 = filter_bp(ch2, 20e6 ,80e6 ,8);
        filtered_signal3 = filter_bp(ch3, 20e6 ,80e6 ,8);

%         filtered_signal1 = preprocess(ch1);
%         filtered_signal2 = preprocess(ch2);
%         filtered_signal3 = preprocess(ch3);

%寻找信号1的所有满足条件的峰值
peaks = find_peaks(filtered_signal1,12);
%遍历所有峰值
for  pi = 1:numel(peaks)
        idx = peaks(pi);
    if idx-(window_length/2-1) <= 0 
        continue;  % 超出范围，执行下一个区间
    end
    if  idx+(window_length/2) > length(filtered_signal1)
        break;  % 超出范围，执行下一个区间
    end
%     取峰值两端一定长度的信号
    signal1 = filtered_signal1(idx-(window_length/2-1):idx+(window_length/2));
    signal2 = filtered_signal2(idx-(window_length/2-1):idx+(window_length/2));
    signal3 = filtered_signal3(idx-(window_length/2-1):idx+(window_length/2));
%     窗口处理
    windows =1:256:length(signal1)-window_length+1;
    for  wi = 1:numel(windows)
        win_signal1 = signal1(windows(wi):windows(wi)+window_length-1);
        win_signal2 = signal2(windows(wi):windows(wi)+window_length-1);
        win_signal3 = signal3(windows(wi):windows(wi)+window_length-1);

%         filtered_signal1 = filter_bp(win_signal1, 20e6 ,80e6 ,8);
%         filtered_signal2 = filter_bp(win_signal2, 20e6 ,80e6 ,8);
%         filtered_signal3 = filter_bp(win_signal3, 20e6 ,80e6 ,8);

%         filtered_signal1 = preprocess(win_signal1);
%         filtered_signal2 = preprocess(win_signal2);
%         filtered_signal3 = preprocess(win_signal3);

        %去直流分量
        signal1_removed = detrend(win_signal1);
        signal2_removed = detrend(win_signal2);
        signal3_removed = detrend(win_signal3);
        % 应用窗函数
        ch1_new = real(windowsignal(signal1_removed));
        ch2_new = real(windowsignal(signal2_removed));
        ch3_new = real(windowsignal(signal3_removed));
        %信号上采样
        ch1_upsp = upsampling(ch1_new,50)';
        ch2_upsp = upsampling(ch2_new,50)';
        ch3_upsp = upsampling(ch3_new,50)';
        %互相关
        [r12_gcc,lags12_gcc] = xcorr(ch1_upsp(:,2),ch2_upsp(:,2),'normalized');
        [r13_gcc,lags13_gcc] = xcorr(ch1_upsp(:,2),ch3_upsp(:,2),'normalized');
        [r23_gcc,lags23_gcc] = xcorr(ch2_upsp(:,2),ch3_upsp(:,2),'normalized');
        R12_gcc = max(r12_gcc);
        R13_gcc = max(r13_gcc);
        R23_gcc = max(r23_gcc);
        t12_gcc = cal_tau(r12_gcc,lags12_gcc');
        t13_gcc = cal_tau(r13_gcc,lags13_gcc');
        t23_gcc = cal_tau(r23_gcc,lags23_gcc');
        Rcorr = (R12_gcc + R13_gcc + R23_gcc) / 3;
        if Rcorr <0.3
            continue;
        end
         %根据互相关曲线的最大值进行数据平移的结果
        ch1_gcc_new = ch1_upsp(:,2);
        ch2_gcc_new = shift_signal(ch2_upsp(:,2),t12_gcc);
        ch3_gcc_new = shift_signal(ch3_upsp(:,2),t13_gcc);
        %设置阈值并寻找峰值
        peaks1= find_max_peaks(ch1_gcc_new,4);
        peaks2= find_max_peaks(ch2_gcc_new,4);
        peaks3= find_max_peaks(ch3_gcc_new,4);
        %匹配峰值得到每一对峰值的x坐标
        matched_peaks_x = match_peaks(peaks1,peaks2,peaks3);
        R12s = []; %用于存储R12s值的空向量
        R13s = []; %用于存储R13s值的空向量
        R23s = []; %用于存储R23s值的空向量

        
        for i = 1 : size(matched_peaks_x, 1)
            %微尺度
            ch1_msw = msw_signal(ch1_gcc_new , matched_peaks_x(i,1) ,400);
            ch2_msw = msw_signal(ch2_gcc_new , matched_peaks_x(i,2) ,400);
            ch3_msw = msw_signal(ch3_gcc_new , matched_peaks_x(i,3) ,400);
            if numel(ch1_msw)~= 800 || numel(ch2_msw)~= 800 || numel(ch3_msw)~= 800
                continue;
            end
            %对微尺度进行互相关
            [R12_msw,lag12_msw] = xcorr(ch1_msw,ch2_msw,'normalized');
            [R13_msw,lag13_msw] = xcorr(ch1_msw,ch3_msw,'normalized');
            [R23_msw,lag23_msw] = xcorr(ch2_msw,ch3_msw,'normalized');
            R12s = [R12s R12_msw];
            R13s = [R13s R13_msw];
            R23s = [R23s R23_msw];
        end
        
        for Ri = 1:size(R12s,2)
%             验证是否三个相关系数都大于0.8
            R12 = max(R12s(: , Ri));
            R13 = max(R13s(: , Ri));
            R23 = max(R23s(: , Ri));

            if R12< 0.5 || R13 < 0.5 || R23 < 0.5
                continue
            end

            t12 = (t12_gcc  + (matched_peaks_x(Ri,1)-matched_peaks_x(Ri,2)))*0.1;
            t13 = (t13_gcc  + (matched_peaks_x(Ri,1)-matched_peaks_x(Ri,3)))*0.1;
            t23 = (t23_gcc  + (matched_peaks_x(Ri,2)-matched_peaks_x(Ri,3)))*0.1;
            
            cos_alpha_0 = c*t23*tand(angle23)/(d*sind(angle23)*(tand(angle23) - tand(angle12))) - c*t12/(d*cosd(angle12)*(tand(angle23)-tand(angle12)));
            cos_beta_0 = (c*t12-d*cos_alpha_0*sind(angle12))/(d*cosd(angle12));
            if abs(cos_beta_0)>1 || abs(cos_alpha_0)>1
               continue;
            end
            x0 = [cos_alpha_0,cos_beta_0];
            % 调用lsqnonlin函数进行优化
            options = optimoptions('lsqnonlin', 'MaxIter', 1000, 'TolFun', 1e-6);
            x = lsqnonlin(@objective, x0, [-1 -1],[1 1], options);
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
           

            % 写入计算后的数据
            fprintf(fileID, '%-13d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
            300000000+idx+windows(wi),matched_peaks_x(Ri,1),matched_peaks_x(Ri,2),matched_peaks_x(Ri,3), t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
        end
    end
end 
% 关闭文件
fclose(fileID);




function delta_t = delta_t(tij,tij_obs)
    delta_t = tij - tij_obs;
end



% 定义计算τij的理想值τ_ij^obs的函数
function tau_ij_obs = calculate_tau_obs(cos_alpha, cos_beta)
    angle12 = -150;
    angle13 = -90;
    angle23 = -30;
    % 使用式(3)计算τij的理想值τ_ij^obs
    tau_ij_obs(1) = (cos_alpha * sind(angle12) + cos_beta * cosd(angle12)) * 20 / 0.299792458;
    tau_ij_obs(2) = (cos_alpha * sind(angle13) + cos_beta * cosd(angle13)) * 20 / 0.299792458;
    tau_ij_obs(3) = (cos_alpha * sind(angle23) + cos_beta * cosd(angle23)) * 20 / 0.299792458;
end


function tau = cal_tau(R, lag)
    % 从数据中找到y的最大值及其索引
    [~, max_index] = max(R);
    tau = lag(max_index,1);
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
        [tau12,R12,lag12] = gccphat(window1,window2, Fs);
        [tau13,R13,lag13] = gccphat(window1,window3, Fs);
        [tau23,R23,lag23] = gccphat(window2,window3, Fs);
        max_R12 = maxvalue(R12);
        max_R13 = maxvalue(R13);
        max_R23 = maxvalue(R23);
        
        % 如果相关系数大于阈值，将窗口添加到结果列表中
        if max_R12 > threshold && max_R13 > threshold && max_R23 > threshold
            correlated_windows = [correlated_windows; i];
        end
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

% 定义目标函数
function F = objective(x)
    % 提取待优化的变量
    cos_alpha = x(1);
    cos_beta = x(2);

    % 计算τij的理想值τ_ij^obs
    tau_ij_obs = calculate_tau_obs(cos_alpha, cos_beta);
    t12 = evalin('base', 't12');
    t13 = evalin('base', 't13');
    t23 = evalin('base', 't23');
    % 计算Δt12, Δt13, Δt23
    delta_t12 = delta_t(t12,tau_ij_obs(1));
    delta_t13 = delta_t(t13,tau_ij_obs(2));
    delta_t23 = delta_t(t23,tau_ij_obs(3));

    % 计算目标函数，即式(4)
    F = (delta_t12^2 + delta_t13^2 + delta_t23^2) / 75;
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


function tau = showfitted(data)
    % 从数据中找到y的最大值及其索引
    [max_value, max_index] = max(data(:, 2));
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
        [max_value_fit, max_index_fit] = max(fit_values_curve);
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


