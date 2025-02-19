function [start_read_loc_chj_top10, top10_r_gccs]= get_match_single_yld_chj(yld_signal_start_loc)
upsampling_factor = 50;
yld_signal_length = 1024;
yld_signal = read_signal('../20240822165932.6610CH1.dat',yld_signal_length,yld_signal_start_loc);

chj_signal_length = 4e4;
chj_signal_start_loc = yld_signal_start_loc+71983665;
chj_signal = read_signal('../2024 822 85933.651462CH1.dat',chj_signal_length,chj_signal_start_loc);

all_locs = [];
all_R_gccs = [];
all_t_gccs = [];
filtered_yld_signal = filter_bp(yld_signal, 20e6 ,80e6 ,5);
processed_yld_signal = real(windowsignal(detrend(filtered_yld_signal)));
yld_signal_up = upsampling(processed_yld_signal, upsampling_factor)';
yld_signal_up = yld_signal_up(:,2);
filtered_chj_signal = filter_bp(chj_signal, 20e6 ,80e6 ,5);
processed_chj_signal = real(windowsignal(detrend(filtered_chj_signal)));
subsignal_length = 4000;
subsignal_start = 1:subsignal_length:length(processed_chj_signal);
for subi = 1:numel(subsignal_start)
    subsignal1 = processed_chj_signal(subsignal_start(subi):subsignal_start(subi)+subsignal_length-1);
    threshold = 0.5 * mean(abs(subsignal1));

    % 寻找峰值
    [peaks, locs] = findpeaks(subsignal1, 'MinPeakHeight', threshold, 'MinPeakDistance', subsignal_length/4);

    % 存储所有峰值和阈值
    all_locs = [all_locs; locs + (subsignal_start(subi) - 1)];
end
% 遍历所有峰值
num_peaks = numel(all_locs);
for pi = 1:num_peaks
    idx = all_locs(pi);
    % 确保峰值不超出信号范围
    if idx - (yld_signal_length / 2 - 1) <= 0 || idx + (yld_signal_length / 2) > length(processed_chj_signal)
        continue;
    end
    % 截取窗口信号
    signal1 = processed_chj_signal(idx-(yld_signal_length/2-1):idx+(yld_signal_length/2));
    % 去直流分量并应用窗函数
    ch1_new = real(windowsignal(detrend(signal1)));
    % 上采样
    ch1_up = upsampling(ch1_new, upsampling_factor)';
    ch1_upsp = ch1_up(:,2);
    %互相关
    [r_gcc,lags_gcc] = xcorr(ch1_upsp,yld_signal_up,'normalized');
    R_gcc = max(r_gcc);
    all_R_gccs = [all_R_gccs;R_gcc];
    t_gcc = cal_tau(r_gcc,lags_gcc');
    all_t_gccs = [all_t_gccs;t_gcc];
end

[top10_r_gccs, top10_r_gcc_idxs] = maxk(all_R_gccs, 10);
start_read_loc_chj_top10 = zeros(1, 10);
for i = 1:10
    start_read_loc_chj_top10(i) = chj_signal_start_loc + all_locs(top10_r_gcc_idxs(i)) - 512 + floor(all_t_gccs(top10_r_gcc_idxs(i)) / 50);
end
end
