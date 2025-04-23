yld_ch1 =read_signal('..\\20240822165932.6610CH1.dat',6e5,469400180);
chj_ch1 =read_signal('..\\20240822165932.6610CH1.dat',6e6,469400180+5000);
chj_ch1 = filter_bp(chj_ch1,30e6,80e6,5);
yld_ch1 = filter_bp(yld_ch1,30e6,80e6,5);



[r_gcc, lags_gcc] = xcorr(yld_ch1, chj_ch1,'none');
t_gcc            = cal_tau(r_gcc, lags_gcc');


[peaks, locs] = findpeaks(chj_ch1, 'MinPeakHeight', 30, 'MinPeakDistance', 256);
all_locs = locs;
% 遍历所有峰值
num_peaks = numel(all_locs);
R_gcc_2 = 0;
T_gcc_2 = 0;
idx_2 = 0;
for pi = 1:num_peaks
    idx = all_locs(pi);
    if idx - (6e5 / 2 - 1) <= 0 || idx + (6e5 / 2) > 6e6
        continue;
    end
    processed_chj_signal1 = chj_ch1(idx - (6e5 / 2)+ 1:idx + (6e5 / 2));
    [r_gcc, lags_gcc] = xcorr(yld_ch1, processed_chj_signal1, 'normalized');
    [R_gcc , max_idx]= max(r_gcc);
    t_gcc_2 = cal_tau(r_gcc, lags_gcc');
    if R_gcc > R_gcc_2
        T_gcc_2 = t_gcc_2;
        idx_2 = idx;
        R_gcc_2 = R_gcc;
    end
end