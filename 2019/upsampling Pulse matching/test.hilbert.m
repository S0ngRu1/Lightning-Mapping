%读取数据
signal_length = 1e8;
ch1 = read_signal('../cross-correlation/20190604164852.7960CH1.dat',signal_length);
ch2 = read_signal('../cross-correlation/20190604164852.7960CH2.dat',signal_length);
ch3 = read_signal('../cross-correlation/20190604164852.7960CH3.dat',signal_length);

fs = 200e6;
N = 3;
d = 20;
c = 0.299792458;
window_length = 1024;

angle12 = -150;
angle13 = -90;
angle23 = -30;

% 打开一个文本文件用于写入运行结果
fileID = fopen('result14.txt', 'w');
% 写入第一行的数据介绍
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
     'Start_loc','Peak1_loc','Peak2_loc','Peak3_loc','t12', 't13', 't23', 'cosα', 'cosβ', 'Azimuth', 'Elevation', 'Rcorr', 't123');




filtered_signal1 = filter_bp(ch1, 20e6 ,80e6 ,8);
filtered_signal2 = filter_bp(ch2, 20e6 ,80e6 ,8);
filtered_signal3 = filter_bp(ch3, 20e6 ,80e6 ,8);

       
%寻找信号1的所有满足条件的峰值
peaks = find_peaks(filtered_signal1,12);
%遍历所有峰值
for  pi = 1:numel(peaks)
        idx = peaks(pi);
    if idx-(window_length/2-1) <= 0 
        continue;  % 超出范围，执行下一个区间
    end
    if  idx+(window_length/2) > length(filtered_signal1)
        break;  % 超出范围，执行下一个区间
    end
%     取峰值两端一定长度的信号
    signal1 = filtered_signal1(idx-(window_length/2-1):idx+(window_length/2));
    signal2 = filtered_signal2(idx-(window_length/2-1):idx+(window_length/2));
    signal3 = filtered_signal3(idx-(window_length/2-1):idx+(window_length/2));
%     窗口处理
    windows =1:256:length(signal1)-window_length+1;
    for  wi = 1:numel(windows)
        win_signal1 = signal1(windows(wi):windows(wi)+window_length-1);
        win_signal2 = signal2(windows(wi):windows(wi)+window_length-1);
        win_signal3 = signal3(windows(wi):windows(wi)+window_length-1);

%         filtered_signal1 = filter_bp(win_signal1, 20e6 ,80e6 ,8);
%         filtered_signal2 = filter_bp(win_signal2, 20e6 ,80e6 ,8);
%         filtered_signal3 = filter_bp(win_signal3, 20e6 ,80e6 ,8);

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
        Rcorr = (R12_gcc + R13_gcc + R23_gcc) / 3;
        if Rcorr <0.3
            continue;
        end
         %根据互相关曲线的最大值进行数据平移的结果
        ch1_gcc_new = ch1_upsp(:,2);
        ch2_gcc_new = shift_signal(ch2_upsp(:,2),t12_gcc);
        ch3_gcc_new = shift_signal(ch3_upsp(:,2),t13_gcc);
        
        %对平移后的信号进行希尔伯特变换
        ch1_ht = imag(hilbert(ch1_gcc_new));
        ch2_ht = imag(hilbert(ch2_gcc_new));
        ch3_ht = imag(hilbert(ch3_gcc_new));

        L1 = length(ch1_ht);
        frequencies1 = (0:L1-1)*(fs/L1); % 计算频率轴
        power_spectrum1 = abs(fft(ch1_ht)).^2 / L1; % 计算功率谱
        L2 = length(ch2_ht);
        frequencies2 = (0:L2-1)*(fs/L2); % 计算频率轴
        power_spectrum2 = abs(fft(ch2_ht)).^2 / L2; % 计算功率谱
        L3 = length(ch3_ht);
        frequencies3 = (0:L3-1)*(fs/L3); % 计算频率轴
        power_spectrum3 = abs(fft(ch3_ht)).^2 / L3; % 计算功率谱


        %设置阈值并寻找峰值
        peaks1= find_max_peaks(power_spectrum1,30);
        peaks2= find_max_peaks(power_spectrum2,30);
        peaks3= find_max_peaks(power_spectrum3,30);
        %匹配峰值得到每一对峰值的x坐标
        matched_peaks_x = match_peaks(peaks1,peaks2,peaks3);
        R12s = []; %用于存储R12s值的空向量
        R13s = []; %用于存储R13s值的空向量
        R23s = []; %用于存储R23s值的空向量

        
        for i = 1 : size(matched_peaks_x, 1)
            %微尺度
            ch1_msw = msw_signal(ch1_gcc_new , matched_peaks_x(i,1) ,400);
            ch2_msw = msw_signal(ch2_gcc_new , matched_peaks_x(i,2) ,400);
            ch3_msw = msw_signal(ch3_gcc_new , matched_peaks_x(i,3) ,400);
            if numel(ch1_msw)~= 800 || numel(ch2_msw)~= 800 || numel(ch3_msw)~= 800
                continue;
            end
            %对微尺度进行互相关
            [R12_msw,lag12_msw] = xcorr(ch1_msw,ch2_msw,'normalized');
            [R13_msw,lag13_msw] = xcorr(ch1_msw,ch3_msw,'normalized');
            [R23_msw,lag23_msw] = xcorr(ch2_msw,ch3_msw,'normalized');
            R12s = [R12s R12_msw];
            R13s = [R13s R13_msw];
            R23s = [R23s R23_msw];
        end
        
        for Ri = 1:size(R12s,2)
%             验证是否三个相关系数都大于0.8
            R12 = max(R12s(: , Ri));
            R13 = max(R13s(: , Ri));
            R23 = max(R23s(: , Ri));

            if R12< 0.5 || R13 < 0.5 || R23 < 0.5
                continue
            end

            t12 = (t12_gcc  + (matched_peaks_x(Ri,1)-matched_peaks_x(Ri,2)))*0.1;
            t13 = (t13_gcc  + (matched_peaks_x(Ri,1)-matched_peaks_x(Ri,3)))*0.1;
            t23 = (t23_gcc  + (matched_peaks_x(Ri,2)-matched_peaks_x(Ri,3)))*0.1;
            
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
            
            t123 = t12 + t23 - t13;
           

            % 写入计算后的数据
            fprintf(fileID, '%-13d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
            400000000+idx+windows(wi),matched_peaks_x(Ri,1),matched_peaks_x(Ri,2),matched_peaks_x(Ri,3), t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
        end
    end
end 
% 关闭文件
fclose(fileID);
