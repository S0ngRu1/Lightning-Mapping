function [start_loc, azimuth, elevation, Rcorr, t123]  = get_2d_result_single_window(start_loc,ch1,ch2,ch3 ,type)
if strcmp(type, 'chj')
    % 从化局
    d12 = 41.6496;
    d13 = 48.5209;
    d23 = 25.0182;
    angle12 = -2.8381;
    angle13 = 28.2006;
    angle23 = 87.3358;
elseif strcmp(type, 'yld')
    % 引雷点
    angle12 = -110.8477;
    angle13 = -65.2405;
    angle23 = -19.6541;
    d12 = 24.9586;
    d13 = 34.9335;
    d23 = 24.9675;

end
start_loc = 0;
azimuth = 0;
elevation = 0;
Rcorr = 0;
t123 = 0;
N = 3;
c = 0.299792458;
fs = 200e6;
window_length = 1024;
bigwindows_length = window_length+100;
window = window_length * 50;
msw_length = 50;
% 在输入信号1上寻峰
% 动态阈值
threshold =  mean(abs(ch1)) + 3*std(ch1);
[peaks, locs] = findpeaks(ch1, 'MinPeakHeight', threshold, 'MinPeakDistance', 256);
% 存储所有峰值
all_peaks = peaks;
all_locs = locs;
all_azimuths = [];
all_t123s = [];
all_elevations = [];
all_Rcorrs = [];
% 遍历所有峰值
num_peaks = numel(all_peaks);
for pi = 1:num_peaks
    idx = all_locs(pi);
    % 确保峰值不超出信号范围
    if idx - (bigwindows_length / 2 - 1) <= 0 || idx + (bigwindows_length / 2) > length(ch1)
        continue;
    end
    [ch1_new, ch2_new, ch3_new] = deal(...
        ch1(idx-(bigwindows_length/2-1):idx+(bigwindows_length/2)), ...
        ch2(idx-(bigwindows_length/2-1):idx+(bigwindows_length/2)), ...
        ch3(idx-(bigwindows_length/2-1):idx+(bigwindows_length/2)));
    [ch1_up, ch2_up, ch3_up] = deal(...
        upsampling(ch1_new, 50)', ...
        upsampling(ch2_new, 50)', ...
        upsampling(ch3_new, 50)');
    ch1_upsp = ch1_up(:,2);
    ch2_upsp = ch2_up(:,2);
    ch3_upsp = ch3_up(:,2);


    ch1_new = ch1_upsp(bigwindows_length*50/2-(window/2-1):bigwindows_length*50/2+(window/2));
    ch2_new = ch2_upsp(bigwindows_length*50/2-(window/2-1):bigwindows_length*50/2+(window/2));
    ch3_new = ch3_upsp(bigwindows_length*50/2-(window/2-1):bigwindows_length*50/2+(window/2));
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

    ismsw = 0;
    if R12_gcc > 0.8 && R13_gcc > 0.8
        %如果R大于某个值则将根据互相关曲线的最大值进行数据平移的结果
        shifted_ch1 = ch1_upsp(bigwindows_length*50/2-(window/2-1):bigwindows_length*50/2+(window/2));
        shifted_ch2 = shift_signal(ch2_upsp,t12_gcc);
        shifted_ch2 = shifted_ch2(bigwindows_length*50/2-(window/2-1):bigwindows_length*50/2+(window/2));
        shifted_ch3 = shift_signal(ch3_upsp,t13_gcc);
        shifted_ch3 = shifted_ch3(bigwindows_length*50/2-(window/2-1):bigwindows_length*50/2+(window/2));

        threshold1 = mean(shifted_ch1);
        threshold2 = mean(shifted_ch2);
        threshold3 = mean(shifted_ch3);
        peaks1= find_peaks(shifted_ch1,threshold1);
        peaks2= find_peaks(shifted_ch2,threshold2);
        peaks3= find_peaks(shifted_ch3,threshold3);

        matched_peaks_x = match_peaks(peaks1,peaks2,peaks3);
        R12s = []; %用于存储R12s值的空向量
        R13s = []; %用于存储R13s值的空向量
        R23s = []; %用于存储R23s值的空向量
        t12s = [];
        t13s = [];
        t23s = [];
        if size(matched_peaks_x,1) ~= 0
            for i = 1 : size(matched_peaks_x, 1)
                %微尺度
                ch1_msw = msw_signal(shifted_ch1 , matched_peaks_x(i,1) ,msw_length,window);
                ch2_msw = msw_signal(shifted_ch2 , matched_peaks_x(i,2) ,msw_length,window);
                ch3_msw = msw_signal(shifted_ch3 , matched_peaks_x(i,3) ,msw_length,window);
                if numel(ch1_msw)~= msw_length*2 || numel(ch2_msw)~= msw_length*2 || numel(ch3_msw)~= msw_length*2
                    continue;
                end
                %对微尺度进行互相关
                [R12_msw,lags12_msw] = xcorr(ch1_msw,ch2_msw,'normalized');
                [R13_msw,lags13_msw] = xcorr(ch1_msw,ch3_msw,'normalized');
                [R23_msw,lags23_msw] = xcorr(ch2_msw,ch3_msw,'normalized');
                if max(R12_msw) > 0.8 && max(R13_msw) > 0.8 && max(R23_msw) > 0.8
                    t12_msw = cal_tau(R12_msw,lags12_msw');
                    t13_msw = cal_tau(R13_msw,lags13_msw');
                    t23_msw = cal_tau(R23_msw,lags23_msw');
                    R12s = [R12s R12_msw];
                    R13s = [R13s R13_msw];
                    R23s = [R23s R23_msw];
                    t12s = [t12s t12_msw];
                    t13s = [t13s t13_msw];
                    t23s = [t23s t23_msw];
                end
            end
            if size(t12s,1)~=0 && size(t13s,1)~=0 && size(t23s,1)~=0
                %                                 %从化局
                t12 = (t12_gcc + mean(t12s))*0.1;
                t13 = (t13_gcc + mean(t13s))*0.1+2;
                t23 = (t23_gcc + mean(t23s))*0.1+2;

                cos_beta_0 =((c*t13*d12*sind(angle12))-(c*t12*sind(angle13)*d13))/(d13*d12*sind(angle12-angle13)) ;
                cos_alpha_0 = ((c*t12)/d12-cos_beta_0*cosd(angle12))/sind(angle12);
                if abs(cos_beta_0)>1 || abs(cos_alpha_0)>1
                    continue;
                end
                x0 = [cos_alpha_0,cos_beta_0];
                % 调用lsqnonlin函数进行优化
                options = optimoptions('lsqnonlin', 'MaxIter', 1000, 'TolFun', 1e-6);
                x = lsqnonlin(@(x) objective(x, t12, t13, t23, type), x0, [-1 -1],[1 1], options);
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
                all_Rcorrs = [all_Rcorrs Rcorr];
                all_azimuths = [all_azimuths Az_deg];
                all_elevations = [all_elevations El_deg];
                all_t123s = [all_t123s t123];
                ismsw = ismsw + 1;
            end
        end
    end
    if ismsw == 0
        %                 从化局
        t12 = t12_gcc *0.1;
        t13 = t13_gcc *0.1+2;
        t23 = t23_gcc *0.1+2;
        cos_beta_0 =((c*t13*d12*sind(angle12))-(c*t12*sind(angle13)*d13))/(d13*d12*sind(angle12-angle13)) ;
        cos_alpha_0 = ((c*t12)/d12-cos_beta_0*cosd(angle12))/sind(angle12);
        if abs(cos_beta_0)>1 || abs(cos_alpha_0)>1
            continue;
        end
        x0 = [cos_alpha_0,cos_beta_0];
        % 调用lsqnonlin函数进行优化
        options = optimoptions('lsqnonlin', 'MaxIter', 1000, 'TolFun', 1e-6);
        x = lsqnonlin(@(x) objective(x, t12, t13, t23,type), x0, [-1 -1],[1 1], options);
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
        all_Rcorrs = [all_Rcorrs Rcorr];
        all_azimuths = [all_azimuths Az_deg];
        all_elevations = [all_elevations El_deg];
        all_t123s = [all_t123s t123];
    end
end

if ~isempty(all_azimuths)
    % 对all_t123 进行升序排序，取最小的并且将对应的取值赋给start_loc, azimuth, elevation, Rcorr, t123
    [~, start_loc] = min(all_t123s);
    [azimuth, elevation, Rcorr, t123] = deal(all_azimuths(start_loc), all_elevations(start_loc), all_Rcorrs(start_loc), all_t123s(start_loc));
    start_loc = 1;
end
end