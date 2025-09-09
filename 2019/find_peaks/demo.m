%读取数据
signal_length = 2e8;
ch1 = read_signal('../cross-correlation/20190604164852.7960CH1.dat',signal_length);
% x = downsample(ch1,50);
% filtered_signalx = preprocess(x);
% plot_signal_spectrum(filtered_signalx);

ch2 = read_signal('../cross-correlation/20190604164852.7960CH2.dat',signal_length);
ch3 = read_signal('../cross-correlation/20190604164852.7960CH3.dat',signal_length);

N = 3;
d = 20;
c = 0.299792458;
window_length = 1024;

angle12 = -150;
angle13 = -90;
angle23 = -30;


% 打开一个文本文件用于写入运行结果
fileID = fopen('result2.preprocess.12.txt', 'w');
% 写入第一行的数据介绍
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
     'Start_loc','Peak_loc','t12', 't13', 't23', 'cosα', 'cosβ', 'Azimuth', 'Elevation', 'Rcorr', 't123');

filtered_signal1 = preprocess(ch1);
filtered_signal2 = preprocess(ch2);
filtered_signal3 = preprocess(ch3);

%寻找信号1的所有满足条件的峰值
peaks = find_peaks(filtered_signal1,12);
%遍历所有峰值
for  pi = 1:numel(peaks)
    idx = peaks(pi);
    if idx-(window_length*1/2-1) <= 0 
        continue;  % 超出范围，执行下一个区间
    end
    if  idx+(window_length*1/2) > length(filtered_signal1)
        break;  % 超出范围，执行下一个区间
    end
%     取峰值两端一定长度的信号
    signal1 = filtered_signal1(idx-(window_length*1/2-1):idx+(window_length*1/2));
    signal2 = filtered_signal2(idx-(window_length*1/2-1):idx+(window_length*1/2));
    signal3 = filtered_signal3(idx-(window_length*1/2-1):idx+(window_length*1/2));
%     窗口处理
    windows =1:256:length(signal1)-window_length+1;
    for  wi = 1:numel(windows)
        win_signal1 = signal1(windows(wi):windows(wi)+window_length-1);
        win_signal2 = signal2(windows(wi):windows(wi)+window_length-1);
        win_signal3 = signal3(windows(wi):windows(wi)+window_length-1);
        
% %         对信号进行滤波处理
%         filtered_signal1 = real(filter_fft(win_signal1, 20e6 ,80e6 ));
%         filtered_signal2 = real(filter_fft(win_signal2, 20e6 ,80e6 ));
%         filtered_signal3 = real(filter_fft(win_signal3, 20e6 ,80e6 ));
%          
%         %带通滤波·
%         filtered_signal1 = filter_bp(win_signal1,20e6,80e6,8);
%         filtered_signal2 = filter_bp(win_signal2,20e6,80e6,8);
%         filtered_signal3 = filter_bp(win_signal3,20e6,80e6,8);

%         filtered_signal1 = preprocess(win_signal1);
%         filtered_signal2 = preprocess(win_signal2);
%         filtered_signal3 = preprocess(win_signal3);

        %去直流分量
        signal1_removed = detrend(win_signal1);
        signal2_removed = detrend(win_signal2);
        signal3_removed = detrend(win_signal3);

        % 对滤波后的信号应用窗函数
        windowed_signal1 = real(windowsignal(signal1_removed));
        windowed_signal2 = real(windowsignal(signal2_removed));
        windowed_signal3 = real(windowsignal(signal3_removed));
        %处理后的信号
        ch1_new = windowed_signal1;
        ch2_new = windowed_signal2;
        ch3_new = windowed_signal3;
        
        % 互相关
        [r12,lags12] = xcorr(ch1_new,ch2_new,'normalized');
        [r13,lags13] = xcorr(ch1_new,ch3_new,'normalized');
        [r23,lags23] = xcorr(ch2_new,ch3_new,'normalized');

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
    
%         % 构建矩阵 A 和向量 B
%         A = [sqrt(3)/2 1/2; sqrt(3)/2 -1/2; 0 1];
%         B = [c*t12/d; c*t13/d; c*t23/d];
%         % 使用左除运算符求解线性方程组的最优解
%         result = A \ B;
%         % 输出最优的cos(α)和cos(β)值
%         cos_alpha_opt = result(1);
%         cos_beta_opt = result(2);
%         if abs(cos_alpha_opt)>1 || abs(cos_beta_opt)>1
%             continue;
%         end
%         Az = atan2( cos_alpha_opt,cos_beta_opt);
%         if abs(cos_beta_opt/cos(Az)) > 1
%             continue;
%         end
%         El = acos( cos_beta_opt/cos(Az) );
%         % 将弧度转换为角度
%         Az_deg = rad2deg(Az);
%         El_deg = rad2deg(El);
%         if Az_deg < 0
%            Az_deg = Az_deg + 360;
%         end


            cos_alpha_0 = c*t23*tand(angle23)/(d*sind(angle23)*(tand(angle23) - tand(angle12))) - c*t12/(d*cosd(angle12)*(tand(angle23)-tand(angle12)));
            cos_beta_0 = (c*t12-d*cos_alpha_0*sind(angle12))/(d*cosd(angle12));
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

        peak_loc = find_max_peaks(ch1_new);
        t123 = t12 + t23 - t13;
        Rcorr = (R12 + R13 + R23)/3;
        
        
        if abs(t123) > 1
           continue;
        end
        % 写入计算后的数据
        fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
             300000000+idx+windows(wi),peak_loc, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
    end
end
% 关闭文件
fclose(fileID);
