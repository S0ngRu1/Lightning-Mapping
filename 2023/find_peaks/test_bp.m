r_loction = 4401e5;
signal_length = 2e5;
ch1 = read_signal('..\\20230718175104.9180CH1.dat',signal_length,r_loction);
ch2 = read_signal('..\\20230718175104.9180CH2.dat',signal_length,r_loction);
ch3 = read_signal('..\\20230718175104.9180CH3.dat',signal_length,r_loction);

filtered_signal1 = filter_bp(ch1, 20e6 ,80e6 ,5);
filtered_signal2 = filter_bp(ch2, 20e6 ,80e6 ,5);
filtered_signal3 = filter_bp(ch3, 20e6 ,80e6 ,5);

window_length = 1024;
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
fileID = fopen('result_text.txt', 'w');
% 写入第一行的数据介绍
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','peak','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');


% 设置动态阈值
subsignal_length = 4000;
subsignal_start = 1:subsignal_length:length(filtered_signal1);
for subi = 1:numel(subsignal_start)
subsignal1 = filtered_signal1(subsignal_start(subi):subsignal_start(subi)+subsignal_length-1);
threshold = 0.5 * mean(abs(subsignal1));
%寻找信号1的所有满足条件的峰值
peaks = find_peaks(filtered_signal1,threshold);

t = 1;
% 遍历所有峰值
for  pi = 1:numel(peaks)
    if peaks(pi) - peaks(t) < 256 && pi ~= 1
        continue;
    end
    idx = peaks(pi);

    if idx-(window_length/2-1) <= 0
        continue;  % 超出范围，执行下一个区间
    end
    if  idx+(window_length/2) > length(filtered_signal1)
        break;  % 超出范围，执行下一个区间
    end
    %取峰值两端一定长度的信号
    signal1 = filtered_signal1(idx-(window_length/2-1):idx+(window_length/2));
    signal2 = filtered_signal2(idx-(window_length/2-1):idx+(window_length/2));
    signal3 = filtered_signal3(idx-(window_length/2-1):idx+(window_length/2));
    

        %去直流分量
        signal1_removed = detrend(signal1);
        signal2_removed = detrend(signal2);
        signal3_removed = detrend(signal3);
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
            r_loction+idx-512,513, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
t = pi;
end
end  
% 关闭文件
fclose(fileID);

