% 参数设置
r_loction = 5e8;
signal_length = 1e8;
windows_length = 1024;
upsampling_factor = 50;
N = 3;
d12 = 24.96;
d13 = 34.93;
d23 = 24.98;
c = 0.299552816;
fs = 200e6;
angle12 = -110.85;
angle13 = -65.24;
angle23 = -19.65;

% 读取信号数据
ch1 = read_signal('..\\20230718175104.9180CH1.dat', signal_length,r_loction);
ch2 = read_signal('..\\20230718175104.9180CH2.dat', signal_length,r_loction);
ch3 = read_signal('..\\20230718175104.9180CH3.dat', signal_length,r_loction);

% 进行带通滤波
filtered_signal1 = filter_bp(ch1, 20e6, 80e6, 5);
filtered_signal2 = filter_bp(ch2, 20e6, 80e6, 5);
filtered_signal3 = filter_bp(ch3, 20e6, 80e6, 5);

% 打开文件用于写入
fileID = fopen('result2_bp_5-6e8_1024_256.txt', 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc', 'peak', 't12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt', 'Azimuth', 'Elevation', 'Rcorr', 't123');

% 计算信号的绝对值均值和标准差
mean_abs_signal1 = mean(abs(filtered_signal1));
std_signal1 = std(filtered_signal1);
% 设置动态阈值
threshold = 0.01*mean_abs_signal1;
% 平滑信号以减少噪声
smoothed_signal1 = movmean(filtered_signal1, 10);
smoothed_signal2 = movmean(filtered_signal2, 10);
smoothed_signal3 = movmean(filtered_signal3, 10);
% 寻找峰值
% 寻找峰值
[peaks, locs] = findpeaks(smoothed_signal1, 'MinPeakHeight', threshold, 'MinPeakDistance', windows_length/4);


% 遍历所有峰值
num_peaks = numel(peaks);
for pi = 1:num_peaks
    idx = locs(pi);

    % 确保峰值不超出信号范围
    if idx - (windows_length / 2 - 1) <= 0 || idx + (windows_length / 2) > length(smoothed_signal1)
        continue;
    end

    % 截取窗口信号
    [signal1, signal2, signal3] = deal(...
        smoothed_signal1(idx-(windows_length/2-1):idx+(windows_length/2)), ...
        smoothed_signal2(idx-(windows_length/2-1):idx+(windows_length/2)), ...
        smoothed_signal3(idx-(windows_length/2-1):idx+(windows_length/2)));


    % 去直流分量并应用窗函数
    [ch1_new, ch2_new, ch3_new] = deal(...
        real(windowsignal(detrend(signal1))), ...
        real(windowsignal(detrend(signal2))), ...
        real(windowsignal(detrend(signal3))));

    % 上采样
      [ch1_upsp, ch2_upsp, ch3_upsp] = deal(...
        upsampling(ch1_new, upsampling_factor)', ...
        upsampling(ch2_new, upsampling_factor)', ...
        upsampling(ch3_new, upsampling_factor)');


    % 计算互相关
    [r12_gcc, lags12_gcc] = xcorr(ch1_upsp(:,2), ch2_upsp(:,2), 'normalized');
    [r13_gcc, lags13_gcc] = xcorr(ch1_upsp(:,2), ch3_upsp(:,2), 'normalized');
    [r23_gcc, lags23_gcc] = xcorr(ch2_upsp(:,2), ch3_upsp(:,2), 'normalized');


    % 计算时间延迟
    [t12, t13, t23] = deal(...
        cal_tau(r12_gcc, lags12_gcc') * 0.1, ...
        cal_tau(r13_gcc, lags13_gcc') * 0.1, ...
        cal_tau(r23_gcc, lags23_gcc') * 0.1);

    % 计算初始角度
    cos_beta_0 = ((c * t13 * d12 * sind(angle12)) - (c * t12 * sind(angle13) * d13)) / (d13 * d12 * sind(angle12 - angle13));
    cos_alpha_0 = ((c * t12) / d12 - cos_beta_0 * sind(angle12)) / sind(angle12);

    % 检查初始角度是否有效
    if abs(cos_beta_0) > 1 || abs(cos_alpha_0) > 1
        continue;
    end

    % 优化
    x0 = [cos_alpha_0, cos_beta_0];
    options = optimoptions('lsqnonlin', 'MaxIter', 1000, 'TolFun', 1e-6);
    x = lsqnonlin(@objective, x0, [-1 -1], [1 1], options);
    [cos_alpha_opt, cos_beta_opt] = deal(x(1), x(2));

    % 检查优化结果
    if abs(cos_alpha_opt) > 1 || abs(cos_beta_opt) > 1
        continue;
    end

    % 计算方位角和俯仰角
    Az = atan2(cos_alpha_opt, cos_beta_opt);
    if abs(cos_beta_opt / cos(Az)) > 1
        continue;
    end
    El = acos(cos_beta_opt / cos(Az));
    [Az_deg, El_deg] = deal(rad2deg(Az), rad2deg(El));
    if Az_deg < 0
        Az_deg = Az_deg + 360;
    end

    % 计算 Rcorr
    R12_gcc = max(xcorr(ch1_upsp(:,2), ch2_upsp(:,2), 'normalized'));
    R13_gcc = max(xcorr(ch1_upsp(:,2), ch3_upsp(:,2), 'normalized'));
    R23_gcc = max(xcorr(ch2_upsp(:,2), ch3_upsp(:,2), 'normalized'));
    Rcorr = (R12_gcc + R13_gcc + R23_gcc) / 3;

    % 计算 t123
    t123 = t12 + t23 - t13;

    % 写入结果
    fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
        r_loction + idx - 512, 513, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr, t123);
end

% 关闭文件
fclose(fileID);

