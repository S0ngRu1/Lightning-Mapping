%读取数据
signal_length = 3000;
ch1 = read_signal('20190604164852.7960CH1.dat',signal_length);
ch2 = read_signal('20190604164852.7960CH2.dat',signal_length);
ch3 = read_signal('20190604164852.7960CH3.dat',signal_length);
plot(ch1)
% 
% x_down = downsample(ch1,20);
% plot(x_down); 
window_length = 1024;
windows =1:256:signal_length-window_length+1;
N = 3;
d = 20;
c = 0.299552816;
fs = 200e6;
% 打开一个文本文件用于写入运行结果
fileID = fopen('result6.txt', 'w');
% 写入第一行的数据介绍
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'corr12', 'corr13', 'corr23');


%遍历所有窗口
for  wi = 1:numel(windows)
%   读信号
    signal1 = ch1(windows(wi):windows(wi)+window_length-1);
    signal2 = ch2(windows(wi):windows(wi)+window_length-1);
    signal3 = ch3(windows(wi):windows(wi)+window_length-1);



%   EEMD分解
%     filtered_eemd1 = filter_eemd(signal1);
%     filtered_eemd2 = filter_eemd(signal2);
%     filtered_eemd3 = filter_eemd(signal3);
%   带通滤波
    filtered_signal1 = bandpass(signal1,[20e6 80e6],200e6);
    filtered_signal2 = bandpass(signal2,[20e6 80e6],200e6);
    filtered_signal3 = bandpass(signal3,[20e6 80e6],200e6);

%     filtered_signal1 = filter_fft(signal1,20e6 ,80e6);
%     filtered_signal2 = filter_fft(signal2,20e6 ,80e6);
%     filtered_signal3 = filter_fft(signal3,20e6 ,80e6);

    %去直流分量
    signal1_removed = detrend(filtered_signal1);
    signal2_removed = detrend(filtered_signal2);
    signal3_removed = detrend(filtered_signal3);

    % 对滤波后的信号应用窗函数
    windowed_signal1 = windowsignal(signal1_removed);
    windowed_signal2 = windowsignal(signal2_removed);
    windowed_signal3 = windowsignal(signal3_removed);
    %处理后的信号
    ch1_new = windowed_signal1;
    ch2_new = windowed_signal2;
    ch3_new = windowed_signal3;
 
    [~,r12,lags12] = gccphat(ch1_new,ch2_new);
    
    [~,r13,lags13] = gccphat(ch1_new,ch3_new);
    [~,r23,lags23] = gccphat(ch2_new,ch3_new);
    
    R12 = max(r12);
    R13 = max(r13);
    R23 = max(r23);
    %对相关系数函数进行上采样
    r12_upsp = upsampling_gc(r12,lags12,8);
    r13_upsp = upsampling_gc(r13,lags13,8);
    r23_upsp = upsampling_gc(r23,lags23,8);

    t12 = showfitted(r12_upsp)*5;
    t13 = showfitted(r13_upsp)*5;
    t23 = showfitted(r23_upsp)*5;

   % 构建矩阵 A 和向量 B
    A = [sqrt(3)/2 1/2; sqrt(3)/2 -1/2; 0 1];
    B = [c*t12/d; c*t13/d; c*t23/d];
    
    % 使用左除运算符求解线性方程组的最优解
    result = A \ B;

    % 输出最优的cos(α)和cos(β)值
    cos_alpha_opt = result(1);
    cos_beta_opt = result(2);
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
    % 写入计算后的数据
        fprintf(fileID, '%-13d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
             400000000+windows(wi), t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, R12,R13,R23);
end
% 关闭文件
fclose(fileID);



