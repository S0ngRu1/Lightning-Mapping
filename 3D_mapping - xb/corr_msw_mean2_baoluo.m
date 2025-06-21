clear;
N = 3;
c = 0.299792458;
fs = 200e6;
upsampling_factor = 50;
window_length = 512;
window = window_length * upsampling_factor;

% 引雷点几何参数
signal_length = 4e7;
r_loction = 3.6e8;
d12 = 24.9586;
d13 = 34.9335;
d23 = 24.9675;
angle12 = -110.8477;
angle13 = -65.2405;
angle23 = -19.6541;

% 读取和滤波数据
ch1 = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction);
ch2 = read_signal('..\\20240822165932.6610CH2.dat',signal_length,r_loction);
ch3 = read_signal('..\\20240822165932.6610CH3.dat',signal_length,r_loction);
filtered_signal1 = filter_bp(ch1,30e6,80e6,5);
filtered_signal2 = filter_bp(ch2,30e6,80e6,5);
filtered_signal3 = filter_bp(ch3,30e6,80e6,5);

% 打开一个文本文件用于写入运行结果
fileID = fopen('result_yld_window512_128_阈值_2n_3.7-3.9_envelope——30——80.txt', 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','peak','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');

% --- 脉冲检测部分 (与您代码相同) ---
% 计算三个通道的包络用于寻峰
env1_full = abs(hilbert(filtered_signal1));
env2_full = abs(hilbert(filtered_signal2));
env3_full = abs(hilbert(filtered_signal3));
% 合成信号
combined_signal = env1_full + env2_full + env3_full;
% 在合成信号上寻峰
noise_level = median(combined_signal) / 0.6745;
dynamic_threshold = 2 * noise_level;
[peaks, locs] = findpeaks(combined_signal, 'MinPeakHeight', dynamic_threshold, 'MinPeakDistance', window_length/4);

all_peaks = peaks;
all_locs = locs;

% 遍历所有峰值
num_peaks = numel(all_peaks);
h = waitbar(0, '正在处理峰值...');
for pi = 1:num_peaks
    waitbar(pi / num_peaks, h, sprintf('正在处理峰值 %d/%d', pi, num_peaks));
    idx = all_locs(pi);
    
    % 确保峰值不超出信号范围
    if idx - (window_length / 2 - 1) <= 0 || idx + (window_length / 2) > length(filtered_signal1)
        continue;
    end
    
    % 截取窗口内的原始射频信号
    [signal1, signal2, signal3] = deal(...
        filtered_signal1(idx-(window_length/2-1):idx+(window_length/2)), ...
        filtered_signal2(idx-(window_length/2-1):idx+(window_length/2)), ...
        filtered_signal3(idx-(window_length/2-1):idx+(window_length/2)));
        
    % ==================================================================
    %                       --- 核心修改处 ---
    % ==================================================================
    % 1. 计算截取出的窗口信号的包络
    env_signal1 = abs(hilbert(signal1));
    env_signal2 = abs(hilbert(signal2));
    env_signal3 = abs(hilbert(signal3));

    % 2. 后续所有处理都基于包络信号进行
    % 去直流分量并应用窗函数 (应用于包络)
    [ch1_new, ch2_new, ch3_new] = deal(...
        real(windowsignal(detrend(env_signal1))), ...
        real(windowsignal(detrend(env_signal2))), ...
        real(windowsignal(detrend(env_signal3))));
    % ==================================================================

    % 上采样 (对处理后的包络进行上采样)
    [ch1_up, ch2_up, ch3_up] = deal(...
        upsampling(ch1_new, upsampling_factor)', ...
        upsampling(ch2_new, upsampling_factor)', ...
        upsampling(ch3_new, upsampling_factor)');
    ch1_upsp = ch1_up(:,2);
    ch2_upsp = ch2_up(:,2);
    ch3_upsp = ch3_up(:,2);
    
    % 互相关 (对上采样后的包络进行互相关)
    [r12_gcc,lags12_gcc] = xcorr(ch1_upsp ,ch2_upsp,'normalized');
    [r13_gcc,lags13_gcc] = xcorr(ch1_upsp ,ch3_upsp,'normalized');
    [r23_gcc,lags23_gcc] = xcorr(ch2_upsp ,ch3_upsp,'normalized');
    
    R12_gcc = max(r12_gcc);
    R13_gcc = max(r13_gcc);
    R23_gcc = max(r23_gcc);
    
    % 假设您已将 cal_tau 升级为高精度版本，例如抛物线插值
    t12_gcc = cal_tau(r12_gcc,lags12_gcc');
    t13_gcc = cal_tau(r13_gcc,lags13_gcc');
    t23_gcc = cal_tau(r23_gcc,lags23_gcc');
    
    % --- 后续的TDOA计算和几何定位部分完全不变 ---
    t12 = t12_gcc *0.1;
    t13 = t13_gcc *0.1;
    t23 = t23_gcc *0.1;
    
    cos_beta_0 =((c*t13*d12*sind(angle12))-(c*t12*sind(angle13)*d13))/(d13*d12*sind(angle12-angle13)) ;
    cos_alpha_0 = ((c*t12)/d12-cos_beta_0*cosd(angle12))/sind(angle12);
    
    if abs(cos_beta_0)>1 || abs(cos_alpha_0)>1
        continue;
    end
    
    x0 = [cos_alpha_0,cos_beta_0];
    options = optimoptions('lsqnonlin', 'Display','off', 'MaxIter', 1000, 'TolFun', 1e-6);
    % 假设您已将 objective 升级为返回向量的正确版本
    x = lsqnonlin(@(x) objective(x, t12, t13, t23,'yld'), x0, [-1 -1],[1 1], options);
    
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
    
    Az_deg = rad2deg(Az);
    El_deg = rad2deg(El);
    if Az_deg < 0
        Az_deg = Az_deg + 360;
    end
    
    t123 = t12 + t23 - t13;
    Rcorr = (R12_gcc + R13_gcc + R23_gcc) / 3;
    
    % 写入计算后的数据
    fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
        r_loction+idx-window/100,window/100, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
end

% 关闭文件和进度条
fclose(fileID);
close(h);