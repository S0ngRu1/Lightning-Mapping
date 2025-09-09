signal_length = 2e5;
r_loction = 4401e5;

ch1 = read_signal('..\\20230718175104.9180CH1.dat',signal_length,r_loction);
ch2 = read_signal('..\\20230718175104.9180CH2.dat',signal_length,r_loction);
ch3 = read_signal('..\\20230718175104.9180CH3.dat',signal_length,r_loction);

filtered_signal1 = filter_bp(ch1, 20e6 ,80e6 ,5);
filtered_signal2 = filter_bp(ch2, 20e6 ,80e6 ,5);
filtered_signal3 = filter_bp(ch3, 20e6 ,80e6 ,5);

%去直流分量
signal1_removed = detrend(filtered_signal1);
signal2_removed = detrend(filtered_signal2);
signal3_removed = detrend(filtered_signal3);

N = 3;
d12 = 24.96;
d13 = 34.93;
d23 = 24.98;
c = 0.299552816;
fs = 200e6;
upsampling_factor = 50;
window_length = 256;
% bigwindows_length = window_length+100;
% window = window_length * upsampling_factor;
% msw_length = 50;

angle12 = -110.85;
angle13 = -65.24;
angle23 = -19.65;

% 打开一个文本文件用于写入运行结果
fileID = fopen('result1_w256_128.txt', 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','peak','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');

step_size = 128;
windows = 1:step_size:(signal_length - window_length + 1);
%遍历所有窗口
for  wi = 1:numel(windows)
    if windows(wi) + window_length - 1 <= signal_length
    signal1 = signal1_removed(windows(wi):windows(wi)+window_length-1);
    signal2 = signal2_removed(windows(wi):windows(wi)+window_length-1);
    signal3 = signal3_removed(windows(wi):windows(wi)+window_length-1);
    end
    %找最大峰值重新取窗口
    max_index = maxindex(signal1);
    %  确保 max_index 足够大以避免负索引
    if max_index+ windows(wi) - (window_length / 2 ) < 1
        % 处理 max_index 太小的情况
        start_index = 1;
    else
        start_index = max_index + windows(wi) - (window_length / 2 );
    end

    % 确保 max_index 足够小以避免超出信号范围
    if max_index + windows(wi) + (window_length / 2 - 1) > length(signal1_removed)
        % 处理 max_index 太大的情况
        end_index = length(signal1_removed);
    else
        end_index = max_index + windows(wi) + (window_length / 2 - 1);
    end

    % 提取有效的窗口数据
    new_signal1 = signal1_removed(start_index:end_index);
    new_signal2 = signal2_removed(start_index:end_index);
    new_signal3 = signal3_removed(start_index:end_index);

%     new_signal1 = signal1_removed(max_index-(window_length/2-1):max_index+(window_length/2-1));
%     new_signal2 = signal2_removed(max_index-(window_length/2-1):max_index+(window_length/2-1));
%     new_signal3 = signal3_removed(max_index-(window_length/2-1):max_index+(window_length/2-1));
    %加窗函数
    window_signal1 = real(windowsignal(new_signal1));
    window_signal2 = real(windowsignal(new_signal2));
    window_signal3 = real(windowsignal(new_signal3));
    %上采样
    ch1_up = upsampling(window_signal1,50)';
    ch2_up = upsampling(window_signal2,50)';
    ch3_up = upsampling(window_signal3,50)';
    ch1_upsp = ch1_up(:,2);
    ch2_upsp = ch2_up(:,2);
    ch3_upsp = ch3_up(:,2);

    %互相关
    [r12_gcc,lags12_gcc] = xcorr(ch1_upsp,ch2_upsp,'normalized');
    [r13_gcc,lags13_gcc] = xcorr(ch1_upsp,ch3_upsp,'normalized');
    [r23_gcc,lags23_gcc] = xcorr(ch2_upsp,ch3_upsp,'normalized');
    R12_gcc = max(r12_gcc);
    R13_gcc = max(r13_gcc);
    R23_gcc = max(r23_gcc);
    t12_gcc = cal_tau(r12_gcc,lags12_gcc');
    t13_gcc = cal_tau(r13_gcc,lags13_gcc');
    t23_gcc = cal_tau(r23_gcc,lags23_gcc');

    t12 = t12_gcc*0.1;
    t13 = t13_gcc*0.1;
    t23 = t23_gcc*0.1;


    %计算
    cos_beta_0 =((c*t13*d12*sind(angle12))-(c*t12*sind(angle13)*d13))/(d13*d12*sind(angle12-angle13)) ;
    cos_alpha_0 = ((c*t12)/d12-cos_beta_0*cosd(angle12))/sind(angle12);
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
    Rcorr = (R12_gcc + R13_gcc + R23_gcc) / 3;
    % 写入计算后的数据
    fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
        r_loction+max_index-window_length/2,window_length/2, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
end
