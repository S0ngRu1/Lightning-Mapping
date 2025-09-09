signal_length = 2e5;
ch1 = read_signal('..\\20230718175104.9180CH1.dat',signal_length);
ch2 = read_signal('..\\20230718175104.9180CH2.dat',signal_length);
ch3 = read_signal('..\\20230718175104.9180CH3.dat',signal_length);

filtered_signal1 = filter_bp(ch1, 20e6 ,80e6 ,5);
filtered_signal2 = filter_bp(ch2, 20e6 ,80e6 ,5);
filtered_signal3 = filter_bp(ch3, 20e6 ,80e6 ,5);

window_length = 1424;
window = 51200;
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
fileID = fopen('result4_pm3.txt', 'w');
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
    %信号上采样
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
    %在平移后的大窗口里面重新取1024窗口
    ch1_newwin = shifted_ch1(35600-(window/2-1):35600+(window/2));
    ch2_newwin = shifted_ch2(35600-(window/2-1):35600+(window/2));
    ch3_newwin = shifted_ch3(35600-(window/2-1):35600+(window/2));

    %设置阈值并寻找峰值
    threshold1 = mean(ch1_newwin);
    threshold2 = mean(ch2_newwin);
    threshold3 = mean(ch3_newwin);
    peaks1= find_max_peaks(ch1_newwin,threshold1);
    peaks2= find_max_peaks(ch2_newwin,threshold2);
    peaks3= find_max_peaks(ch3_newwin,threshold3);
    %匹配峰值得到每一对峰值的x坐标
    matched_peaks_x = match_peaks(peaks1,peaks2,peaks3);

    R12s = []; %用于存储R12s值的空向量
    R13s = []; %用于存储R13s值的空向量
    R23s = []; %用于存储R23s值的空向量

    tmsw = 1;
    for i = 1 : size(matched_peaks_x, 1)
        if matched_peaks_x(i,1) - matched_peaks_x(tmsw,1) < 2560 && i ~= 1   %寻峰间隔
            continue;
        end
        %微尺度
        ch1_msw = msw_signal(ch1_newwin , matched_peaks_x(i,1) ,200);
        ch2_msw = msw_signal(ch2_newwin , matched_peaks_x(i,2) ,200);
        ch3_msw = msw_signal(ch3_newwin , matched_peaks_x(i,3) ,200);
        if numel(ch1_msw)~= 400 || numel(ch2_msw)~= 400 || numel(ch3_msw)~= 400
            continue;
        end
        %对微尺度进行互相关
        [R12_msw,lag12_msw] = xcorr(ch1_msw,ch2_msw,'normalized');
        [R13_msw,lag13_msw] = xcorr(ch1_msw,ch3_msw,'normalized');
        [R23_msw,lag23_msw] = xcorr(ch2_msw,ch3_msw,'normalized');
        R12s = [R12s R12_msw];
        R13s = [R13s R13_msw];
        R23s = [R23s R23_msw];
        tmsw = i;
    end


    for Ri = 1:size(R12s,2)
        %验证是否三个相关系数都大于0.
        R12 = max(R12s(: , Ri));
        R13 = max(R13s(: , Ri));
        R23 = max(R23s(: , Ri));

        if R12< 0.6 || R13 < 0.6 || R23 < 0.6
            continue
        end
        fitted_peak1  = fitpeak(ch1_newwin,matched_peaks_x(i,1));
        fitted_peak2  = fitpeak(ch2_newwin,matched_peaks_x(i,2));
        fitted_peak3  = fitpeak(ch3_newwin,matched_peaks_x(i,3));
        t12 = (t12_gcc + (fitted_peak1-fitted_peak2))*0.1;
        t13 = (t13_gcc + (fitted_peak1-fitted_peak3))*0.1;
        t23 = (t23_gcc + (fitted_peak2-fitted_peak3))*0.1;

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
            401000000+idx-512,513, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
        t = pi;
    end
end
% 关闭文件
fclose(fileID);

