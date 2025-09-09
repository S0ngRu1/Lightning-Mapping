%读取数据
signal_length = 2e5;
ch1 = read_signal('..\\20230718175104.9180CH1.dat',signal_length);
ch2 = read_signal('..\\20230718175104.9180CH2.dat',signal_length);
ch3 = read_signal('..\\20230718175104.9180CH3.dat',signal_length);


% plot_signal_spectrum(ch1);
% x_down = downsample(ch1,20);
% plot(x_down); 
window_length = 1024;   
windows =1:256:signal_length-window_length+1;
N = 3;
d12 = 24.96;
d13 = 34.93;
d23 = 24.98;
c = 0.299552816;
fs = 200e6;

angle12 = -110.85;
angle13 = -65.24;
angle23 = -19.65;

% 打开一个文本文件用于写入运行结果
fileID = fopen('result1.txt', 'w');
% 写入第一行的数据介绍
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');

%   滤波
    filtered_signal1 = filter_bp(ch1, 20e6 ,80e6 ,8);
    filtered_signal2 = filter_bp(ch2, 20e6 ,80e6 ,8);
    filtered_signal3 = filter_bp(ch3, 20e6 ,80e6 ,8);

%遍历所有窗口
for  wi = 1:numel(windows)
    signal1 = filtered_signal1(windows(wi):windows(wi)+window_length-1);
    signal2 = filtered_signal2(windows(wi):windows(wi)+window_length-1);
    signal3 = filtered_signal3(windows(wi):windows(wi)+window_length-1);

    %去直流分量
    signal1_removed = detrend(signal1);
    signal2_removed = detrend(signal2);
    signal3_removed = detrend(signal3);

    % 应用窗函数
    windowed_signal1 = windowsignal(signal1_removed);
    windowed_signal2 = windowsignal(signal2_removed);
    windowed_signal3 = windowsignal(signal3_removed);

    %处理后的信号
    ch1_new = real(windowed_signal1);
    ch2_new = real(windowed_signal2);
    ch3_new = real(windowed_signal3);

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

    t12 = t12_gcc*0.1;
    t13 = t13_gcc*0.1;
    t23 = t23_gcc*0.1;
%     %对相关系数函数进行上采样
%     r12_upsp = upsampling_gc(r12,lags12,10);
%     r13_upsp = upsampling_gc(r13,lags13,10);
%     r23_upsp = upsampling_gc(r23,lags23,10);
% 
%     t12 = showfitted(r12_upsp)*5;
%     t13 = showfitted(r13_upsp)*5;
%     t23 = showfitted(r23_upsp)*5;

%     % 构建矩阵 A 和向量 B
%     A = [sqrt(3)/2 1/2; sqrt(3)/2 -1/2; 0 1];
%     B = [c*t12/d; c*t13/d; c*t23/d];
%     % 使用左除运算符求解线性方程组的最优解
%     result = A \ B;
%     % 输出最优的cos(α)和cos(β)值
%     cos_alpha_opt = result(1);
%     cos_beta_opt = result(2);

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
        fprintf(fileID, '%-13d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
             440100000+windows(wi), t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
end
% 关闭文件
fclose(fileID);