signal_length = 2e5;
ch1 = read_signal('..\\20230718175104.9180CH1.dat',signal_length);
ch2 = read_signal('..\\20230718175104.9180CH2.dat',signal_length);
ch3 = read_signal('..\\20230718175104.9180CH3.dat',signal_length);

% subplot(3,1,1);plot(ch1);title('ch1');xlabel('采样点数');ylabel('幅值');
% subplot(3,1,2);plot(ch2);title('ch2');xlabel('采样点数');ylabel('幅值');
% subplot(3,1,3);plot(ch3);title('ch3');xlabel('采样点数');ylabel('幅值');

filtered_signal1 = filter_bp(ch1, 20e6 ,80e6 ,5);
filtered_signal2 = filter_bp(ch2, 20e6 ,80e6 ,5);
filtered_signal3 = filter_bp(ch3, 20e6 ,80e6 ,5);

% x = downsample(ch1,50);
% plot(filtered_signal1);


N = 3;
d12 = 24.96;
d13 = 34.93;
d23 = 24.98;
c = 0.299552816;
fs = 200e6;
window_length = 1424;
window = 51200;

angle12 = -110.85;
angle13 = -65.24;
angle23 = -19.65;

% 打开一个文本文件用于写入运行结果
fileID = fopen('result8_pm.txt', 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','peak','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');
%寻找信号1的所有满足条件的峰值
threshold = 0.5 * mean(filtered_signal1);
peaks = find_peaks(filtered_signal1,threshold);

t = 1;
%遍历所有峰值
for  pi = 1:numel(peaks)
    if peaks(pi) - peaks(t) < 256 && pi ~= 1   %寻峰间隔
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
        ch1_win = real(windowsignal(signal1_removed));
        ch2_win = real(windowsignal(signal2_removed));
        ch3_win = real(windowsignal(signal3_removed));

        %上采样
        ch1_up = upsampling(ch1_win,50)';
        ch2_up = upsampling(ch2_win,50)';
        ch3_up = upsampling(ch3_win,50)';
        ch1_upsp = ch1_up(:,2);
        ch2_upsp = ch2_up(:,2);
        ch3_upsp = ch3_up(:,2);

        %取窗口做互相关
        ch1_new = ch1_upsp(35600-(window/2-1):35600+(window/2));
        ch2_new = ch2_upsp(35600-(window/2-1):35600+(window/2));
        ch3_new = ch3_upsp(35600-(window/2-1):35600+(window/2));
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

        %根据互相关曲线的最大值进行数据平移的结果
        shifted_ch1 =  ch1_upsp;
        shifted_ch2 = shift_signal(ch2_upsp,t12_gcc);
        shifted_ch3 = shift_signal(ch3_upsp,t13_gcc);


        subplot(3,1,1);plot(shifted_ch1);title('ch1');xlabel('采样点数');ylabel('幅值');
        subplot(3,1,2);plot(shifted_ch2);title('ch2');xlabel('采样点数');ylabel('幅值');
        subplot(3,1,3);plot(shifted_ch3);title('ch3');xlabel('采样点数');ylabel('幅值');

        max_index = maxindex(shifted_ch1);
        ch1_msw = msw_signal(shifted_ch1 , max_index ,400);
        ch2_msw = msw_signal(shifted_ch2 , max_index ,400);
        ch3_msw = msw_signal(shifted_ch3 , max_index ,400);
        if numel(ch1_msw)~= 800 || numel(ch2_msw)~= 800 || numel(ch3_msw)~= 800
            continue;
        end
        [R12_msw,lag12_msw] = xcorr(ch1_msw,ch2_msw,'normalized');
        [R13_msw,lag13_msw] = xcorr(ch1_msw,ch3_msw,'normalized');
        [R23_msw,lag23_msw] = xcorr(ch2_msw,ch3_msw,'normalized');
        t12_msw = cal_tau(R12_msw,lag12_msw');
        t13_msw = cal_tau(R13_msw,lag13_msw');
        t23_msw = cal_tau(R23_msw,lag23_msw');
        t12 = t12_gcc *0.1 + t12_msw*0.1;
        t13 = t13_gcc *0.1 + t13_msw*0.1;
        t23 = t23_gcc *0.1 + t23_msw*0.1;
        

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
            440100000+idx-512,513, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
t = pi;
end
% 关闭文件
fclose(fileID);

