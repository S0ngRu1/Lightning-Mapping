signal_length = 2e5;
ch1 = read_signal('..\\20230718175104.9180CH1.dat',signal_length);
ch2 = read_signal('..\\20230718175104.9180CH2.dat',signal_length);
ch3 = read_signal('..\\20230718175104.9180CH3.dat',signal_length);

data1 = load('noisy1.mat');
Noise1 = (data1.ch1)';
data2 = load('noisy2.mat');
Noise2 = (data2.ch2)';
data3 = load('noisy3.mat');
Noise3 = (data3.ch3)';
% 测量误差协方差
R1 = cov(Noise1); 
R2 = cov(Noise2);
R3 = cov(Noise3);
% 预测误差比较大的时候
Q = 0.5;

filtered_signal1 = KalmanFilter(ch1,Q,R1);
filtered_signal2 = KalmanFilter(ch2,Q,R2);
filtered_signal3 = KalmanFilter(ch3,Q,R3);

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
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','peak','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');
%寻找信号1的所有满足条件的峰值
% peaks = find_peaks(filtered_signal1,92);
% t = 1;
% %遍历所有峰值
% for  pi = 1:numel(peaks)
%     if peaks(pi) - peaks(t) < 256 && pi ~= 1
%         continue;
%     end
%     idx = peaks(pi);
%     if idx-(window_length/2-1) <= 0
%         continue;  % 超出范围，执行下一个区间
%     end
%     if  idx+(window_length/2) > length(filtered_signal1)
%         break;  % 超出范围，执行下一个区间
%     end

%     %取峰值两端一定长度的信号
%     signal1 = filtered_signal1(idx-(window_length/2-1):idx+(window_length/2));
%     signal2 = filtered_signal2(idx-(window_length/2-1):idx+(window_length/2));
%     signal3 = filtered_signal3(idx-(window_length/2-1):idx+(window_length/2));

    % 遍历所有窗口
    for  wi = 1:numel(windows)
        win_signal1 = filtered_signal1(windows(wi):windows(wi)+window_length-1);
        win_signal2 = filtered_signal2(windows(wi):windows(wi)+window_length-1);
        win_signal3 = filtered_signal3(windows(wi):windows(wi)+window_length-1);

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

        t12 = t12_gcc*0.1;
        t13 = t13_gcc*0.1;
        t23 = t23_gcc*0.1;

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
            440100000+windows(wi),513, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
    end
    t = pi;
% 关闭文件
fclose(fileID);

function [Output, NumSegments] = KFrame(Input, WindowLength, Window, HoppingSize)
% Chopper windows the signal based on window length, shift percantage and
% uses Hamming windowing technique.

% Number of samples to hop.
HoppingSamples = fix(WindowLength.*HoppingSize);

% Number of segments.
NumSegments = fix(((length(Input)-WindowLength)/HoppingSamples) + 1);

% Index matrix which guides the signal through chopping process.
Index = (repmat(1:WindowLength,NumSegments,1) + repmat((0:(NumSegments-1))'*HoppingSamples,1,WindowLength))';

% Final window which multiplies with original signal to give pieces of it.
FinalWindow = repmat(Window,1,NumSegments);

% Ta-da... 
Output = Input(Index).*FinalWindow;
end


function T = KalmanFilter(s,Q,R)
% --------------------------------
% T = KalmanFilter(s,Q,R)
% s 位移，列向量
% Q 状态误差协方差矩
% R 观测误差协方差矩
% --------------------------------
% 建立质点运动方程:
% S(k+1) = 1*S(k) + T*V(k) + 0.5*T^2*a
% V(k+1) = 0*S(k) + 1*V(k) + T*a
% 建立观测方程:
% y(k+1) = S(k+1) + v(k+1)
% 即：
% X(k+1) = A * X(k) + G*w(k+1); 预测模型
% y(k+1) = H * X(k+1) + v(k+1); 观测模型
N = length(s);
T = 1; % 采样间隔默认为1
A = [1 T;0 1]; % 状态转移矩阵
G = [T^2/2;T]; % 控制量矩阵
H = [1 0]; % 观测矩阵
% 初始化第一个状态
Xu = [s(1); 0];
Pu = [0 0;0 0];
I = [1 0;0 1];
T = zeros(N,1);
for i = 2:N
Xp = A * Xu;
Pp = A * Pu * A' + G * Q * G';
K = Pp * H' * ( H * Pp * H' + R)^-1;
Xu = ( I - K * H ) * Xp + K * s(i);
Pu = ( I - K * H ) * Pp;
T(i) = Xu(1);
end
end