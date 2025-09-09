signal_length = 1e5;
r_loction = 4e8;

ch1 = read_signal('..\\20230718175104.9180CH1.dat',signal_length,r_loction);
ch2 = read_signal('..\\20230718175104.9180CH2.dat',signal_length,r_loction);
ch3 = read_signal('..\\20230718175104.9180CH3.dat',signal_length,r_loction);

filtered_signal1 = filter_bp(ch1, 20e6 ,80e6 ,5);
filtered_signal2 = filter_bp(ch2, 20e6 ,80e6 ,5);
filtered_signal3 = filter_bp(ch3, 20e6 ,80e6 ,5);

N = 3;
d12 = 24.96;
d13 = 34.93;
d23 = 24.98;
c = 0.299552816;
fs = 200e6;
upsampling_factor = 50;
window_length = 1024;
bigwindows_length = window_length+100;
window = window_length * upsampling_factor;
msw_length = 50;
onnxPath = 'lightning_cnn.onnx';
angle12 = -110.85;
angle13 = -65.24;
angle23 = -19.65;
% onnxModel = importONNXNetwork(onnxPath, 'OutputLayerType', 'classification');
% 打开一个文本文件用于写入运行结果
fileID = fopen('result_cnn.txt', 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');

% 平滑信号以减少噪声
smoothed_signal1 = movmean(filtered_signal1, 10);
smoothed_signal2 = movmean(filtered_signal2, 10);
smoothed_signal3 = movmean(filtered_signal3, 10);

% 遍历所有信号窗口，使用模型预测是否含有信号
windows =1:256:signal_length-window_length+1;
%遍历所有窗口
for  wi = 1:numel(windows)
    [signal1, signal2, signal3] = deal(...
        smoothed_signal1(windows(wi):windows(wi)+(window_length)-1), ...
        smoothed_signal2(windows(wi):windows(wi)+(window_length)-1), ...
        smoothed_signal3(windows(wi):windows(wi)+(window_length)-1));

    predictedLabel = signal_to_onnx(signal1, onnxPath, fs);
    if predictedLabel == 2

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
        [r12_gcc,lags12_gcc] = xcorr(ch1_upsp,ch2_upsp,'normalized');
        [r13_gcc,lags13_gcc] = xcorr(ch1_upsp,ch3_upsp,'normalized');
        [r23_gcc,lags23_gcc] = xcorr(ch2_upsp,ch3_upsp,'normalized');
        R12_gcc = max(r12_gcc);
        R13_gcc = max(r13_gcc);
        R23_gcc = max(r23_gcc);
        t12_gcc = cal_tau(r12_gcc,lags12_gcc');
        t13_gcc = cal_tau(r13_gcc,lags13_gcc');
        t23_gcc = cal_tau(r23_gcc,lags23_gcc');
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
            r_loction+windows, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
    end
end
% 关闭文件
fclose(fileID);

