signal_length = 2e5;
ch1 = read_signal('..\\20230718175104.9180CH1.dat',signal_length);
ch2 = read_signal('..\\20230718175104.9180CH2.dat',signal_length);
ch3 = read_signal('..\\20230718175104.9180CH3.dat',signal_length);

% filtered_signal1 = filter_bp(ch1, 20e6 ,80e6 ,8);
% filtered_signal2 = filter_bp(ch2, 20e6 ,80e6 ,8);
% filtered_signal3 = filter_bp(ch3, 20e6 ,80e6 ,8);

window_length = 1024;
N = 3;
d12 = 24.96;
d13 = 34.93;
d23 = 24.98;
c = 0.299552816;
fs = 200e6;
filtered_signal1 = bandpass(ch1,[20e6 80e6],fs);
filtered_signal2 = bandpass(ch2,[20e6 80e6],fs);
filtered_signal3 = bandpass(ch3,[20e6 80e6],fs);
angle12 = -110.85;
angle13 = -65.24;
angle23 = -19.65;

% 打开一个文本文件用于写入运行结果
fileID = fopen('result1_envelopxocrr_bp.txt', 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','peak','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');
%寻找信号1的所有满足条件的峰值
peaks = find_peaks(filtered_signal1,4);
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
    signal1 = filtered_signal1(idx-(window_length/2-1):idx+(window_length/2));
    signal2 = filtered_signal2(idx-(window_length/2-1):idx+(window_length/2));
    signal3 = filtered_signal3(idx-(window_length/2-1):idx+(window_length/2));
    windows =1:256:length(signal1)-window_length+1;
    for  wi = 1:numel(windows)
        win_signal1 = signal1(windows(wi):windows(wi)+window_length-1);
        win_signal2 = signal2(windows(wi):windows(wi)+window_length-1);
        win_signal3 = signal3(windows(wi):windows(wi)+window_length-1);
        ch1_new = real(preprocessing(win_signal1));
        ch2_new = real(preprocessing(win_signal2));
        ch3_new = real(preprocessing(win_signal3));
        ch1_upsp = upsampling(ch1_new,8)';
        ch2_upsp = upsampling(ch2_new,8)';
        ch3_upsp = upsampling(ch3_new,8)';
        [r12_gcc,lags12_gcc] = xcorr(ch1_upsp(:,2),ch2_upsp(:,2),'normalized');
        [r13_gcc,lags13_gcc] = xcorr(ch1_upsp(:,2),ch3_upsp(:,2),'normalized');
        [r23_gcc,lags23_gcc] = xcorr(ch2_upsp(:,2),ch3_upsp(:,2),'normalized');
%         [up1, ~] = envelope(ch1_new,30,"peak");
%         [up2, ~] = envelope(ch2_new,30,"peak");
%         [up3, ~] = envelope(ch3_new,30,"peak");
%         [r12_gcc,lags12_gcc] = xcorr(up1,up2,'normalized');
%         [r13_gcc,lags13_gcc] = xcorr(up1,up3,'normalized');
%         [r23_gcc,lags23_gcc] = xcorr(up2,up3,'normalized');


        R12_gcc = max(r12_gcc);
        R13_gcc = max(r13_gcc);
        R23_gcc = max(r23_gcc);
        r12_gcc_up = upsampling_gc(r12_gcc,lags12_gcc,10);
        r13_gcc_up = upsampling_gc(r13_gcc,lags13_gcc,10);
        r23_gcc_up = upsampling_gc(r23_gcc,lags23_gcc,10);
        t12_gcc = showfitted(r12_gcc_up);
        t13_gcc = showfitted(r13_gcc_up);
        t23_gcc = showfitted(r23_gcc_up);
%         t12_gcc = cal_tau(r12_gcc,lags12_gcc');
%         t13_gcc = cal_tau(r13_gcc,lags13_gcc');
%         t23_gcc = cal_tau(r23_gcc,lags23_gcc');
%         ch1_gcc_new = ch1_upsp(:,2);
%         ch2_gcc_new = shift_signal(ch2_upsp(:,2),t12_gcc);
%         ch3_gcc_new = shift_signal(ch3_upsp(:,2),t13_gcc);
%         peaks1= find_peaks(ch1_gcc_new,2);
%         peaks2= find_peaks(ch2_gcc_new,2);
%         peaks3= find_peaks(ch3_gcc_new,2);
%         %匹配峰值得到每一对峰值的x坐标
%         matched_peaks_x = match_peaks(peaks1,peaks2,peaks3);
% 
%         t_msw = 1;
%         for i = 1 : size(matched_peaks_x, 1)
%             if matched_peaks_x(i,1) - matched_peaks_x(t_msw,1) < 60 && t_msw ~= 1   %微尺度窗口间隔
%                 continue;
%             end
% 
%             ch1_msw = msw_signal(ch1_gcc_new , matched_peaks_x(i,1) ,30);
%             ch2_msw = msw_signal(ch2_gcc_new , matched_peaks_x(i,2) ,30);
%             ch3_msw = msw_signal(ch3_gcc_new , matched_peaks_x(i,3) ,30);
%             if numel(ch1_msw)~= 60 || numel(ch2_msw)~= 60 || numel(ch3_msw)~= 60
%                 continue;
%             end
%             [up1, ~] = envelope(ch1_msw,30,"peak");
%             [up2, ~] = envelope(ch2_msw,30,"peak");
%             [up3, ~] = envelope(ch3_msw,30,"peak");
%             peak1_index = find_max_peaks(up1);
%             peak2_index = find_max_peaks(up2);
%             peak3_index = find_max_peaks(up3);
%             if peak1_index == -1 || peak2_index == -1 ||peak3_index == -1
%                 continue;
%             end
%             fitted_peak1  = fitpeak(up1,peak1_index);
%             fitted_peak2  = fitpeak(up2,peak2_index);
%             fitted_peak3  = fitpeak(up3,peak3_index);
%             t12 = (t12_gcc + (fitted_peak1-fitted_peak2))*0.1;
%             t13 = (t13_gcc + (fitted_peak1-fitted_peak3))*0.1;
%             t23 = (t23_gcc + (fitted_peak2-fitted_peak3))*0.1;
            t12 = (t12_gcc)*5/8;
            t13 = (t13_gcc)*5/8;
            t23 = (t23_gcc)*5/8;
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
                440100000+idx+windows(wi),513, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
%             t_msw = i;
%         end
    end
    t = pi;
end
% 关闭文件
fclose(fileID);

