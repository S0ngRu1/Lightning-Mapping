N = 3;
c = 0.299792458;
fs = 200e6;
upsampling_factor = 50;
window_length = 1024;
bigwindows_length = window_length+100;
window = window_length * upsampling_factor;
msw_length = 50;
% % 从化局
% d12 = 41.6496;
% d13 = 48.5209;
% d23 = 25.0182;
% angle12 = -2.8381;
% angle13 = 28.2006;
% angle23 = 87.3358;

% signal_length = 3e8;
% r_loction = 3e8;
% ch1 = read_signal('..\\2024 822 85933.651462CH1.dat',signal_length,r_loction);
% ch2 = read_signal('..\\2024 822 85933.651462CH2.dat',signal_length,r_loction);
% ch3 = read_signal('..\\2024 822 85933.651462CH3.dat',signal_length,r_loction165/5);

% %引雷点
% signal_length = 1e8;
% r_loction = 4.5e8;
d12 = 24.9586;
d13 = 34.9335;
d23 = 24.9675;
angle12 = -110.8477;
angle13 = -65.2405;
angle23 = -19.6541;
% ch1 = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction);
% ch2 = read_signal('..\\20240822165932.6610CH2.dat',signal_length,r_loction);
% ch3 = read_signal('..\\20240822165932.6610CH3.dat',signal_length,r_loction);

%滤波
original_signal_length = 3e8;
original_signal_loc = 3e8;
sub_filter_signal_length = 60000;
yld_signal1 = read_signal('..\\20240822165932.6610CH1.dat', original_signal_length, original_signal_loc);
yld_signal2 = read_signal('..\\20240822165932.6610CH2.dat', original_signal_length, original_signal_loc);
yld_signal3 = read_signal('..\\20240822165932.6610CH3.dat', original_signal_length, original_signal_loc);
filtered_yld_signal1 = rfi_filter(yld_signal1,sub_filter_signal_length);
filtered_yld_signal2 = rfi_filter(yld_signal2,sub_filter_signal_length);
filtered_yld_signal3 = rfi_filter(yld_signal3,sub_filter_signal_length);

% filtered_signal1 = filter_bp(ch1, 30e6 ,80e6 ,5);
% filtered_signal2 = filter_bp(ch2, 30e6 ,80e6 ,5);
% filtered_signal3 = filter_bp(ch3, 30e6 ,80e6 ,5);


noise = read_signal('..\\20240822165932.6610CH1.dat',60000,2e8);
filtered_noise = rfi_filter(noise,60000);
threshold = 7*std(filtered_noise);

% 打开一个文本文件用于写入运行结果
fileID = fopen('result_yld_th1_3-6e8_RFI.txt', 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','peak','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');

% 设置动态阈值
% all_thresholds = [];
% % 设置动态阈值，每4000点取一个阈值
% subsignal_length = 4000;
% subsignal_start = 1:subsignal_length:length(filtered_signal1);
% for subi = 1:numel(subsignal_start)
% subsignal1 = filtered_signal1(subsignal_start(subi):subsignal_start(subi)+subsignal_length-1);
% % threshold = 0.5 * mean(abs(subsignal1));
% threshold =  mean(abs(subsignal1)) + 3*std(subsignal1);

% 寻找峰值
[peaks, locs] = findpeaks(filtered_yld_signal1, 'MinPeakHeight', threshold, 'MinPeakDistance', window_length/4);

% 存储所有峰值和阈值
all_peaks = peaks;
%     all_thresholds = [all_thresholds; threshold];
%     all_locs = [all_locs; locs + (subsignal_start(subi) - 1)];
all_locs = locs;

% 遍历所有峰值
num_peaks = numel(all_peaks);
% 创建进度条
h = waitbar(0, '正在处理峰值...');
for pi = 1:num_peaks
    waitbar(pi / num_peaks, h, sprintf('正在处理峰值 %d/%d', pi, num_peaks));
    idx = all_locs(pi);

    % 确保峰值不超出信号范围
    if idx - (bigwindows_length / 2 - 1) <= 0 || idx + (bigwindows_length / 2) > length(filtered_yld_signal1)
        continue;
    end

    % 截取窗口信号
    [signal1, signal2, signal3] = deal(...
        filtered_yld_signal1(idx-(bigwindows_length/2-1):idx+(bigwindows_length/2)), ...
        filtered_yld_signal2(idx-(bigwindows_length/2-1):idx+(bigwindows_length/2)), ...
        filtered_yld_signal3(idx-(bigwindows_length/2-1):idx+(bigwindows_length/2)));
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
%                                 %从化局
%                                 t12 = (t12_gcc + mean(t12s))*0.1;
%                                 t13 = (t13_gcc + mean(t13s))*0.1+2;
%                                 t23 = (t23_gcc + mean(t23s))*0.1+2;
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
                    original_signal_loc+15000+idx-window/100,window/100, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
                ismsw = ismsw + 1;
            end
        end
    end
    if ismsw == 0
%                 从化局
%                 t12 = t12_gcc *0.1;
%                 t13 = t13_gcc *0.1+2;
%                 t23 = t23_gcc *0.1+2;
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
            original_signal_loc+15000+idx-window/100,window/100, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
    end
end
% 关闭文件
fclose(fileID);
% 关闭进度条
close(h);


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