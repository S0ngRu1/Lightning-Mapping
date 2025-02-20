function [start_read_loc_chj, r_gccs]= get_match_single_yld_chj_siddle(yld_signal_start_loc)
yld_signal_length = 1024;
yld_signal = read_signal('../20240822165932.6610CH1.dat',yld_signal_length,yld_signal_start_loc);
chj_signal_length = 1e8;
chj_signal_start_loc = 4.5e8;
chj_signal = read_signal('../2024 822 85933.651462CH1.dat',chj_signal_length,chj_signal_start_loc);
all_R_gccs = [];
all_t_gccs = [];
filtered_yld_signal = filter_bp(yld_signal, 20e6 ,80e6 ,5);
processed_yld_signal = real(windowsignal(detrend(filtered_yld_signal)));
subsignal_length = yld_signal_length/4;
subsignal_start = 1:subsignal_length:chj_signal_length;
for subi = 1:numel(subsignal_start)
    if subsignal_start(subi)+yld_signal_length-1 > chj_signal_length
        continue
    end
    subsignal1 = chj_signal(subsignal_start(subi):subsignal_start(subi)+yld_signal_length-1);
    filtered_chj_signal = filter_bp(subsignal1, 20e6 ,80e6 ,5);
    processed_chj_signal = real(windowsignal(detrend(filtered_chj_signal)));
    [r_gcc,lags_gcc] = xcorr(processed_chj_signal,processed_yld_signal,'normalized');
    R_gcc = max(r_gcc);
    all_R_gccs = [all_R_gccs;R_gcc];
    t_gcc = cal_tau(r_gcc,lags_gcc');
    all_t_gccs = [all_t_gccs;t_gcc];
end

[r_gccs, r_gcc_idxs] = max(all_R_gccs);
start_read_loc_chj = chj_signal_start_loc + subsignal_start(r_gcc_idxs) + floor(all_t_gccs(r_gcc_idxs));

end