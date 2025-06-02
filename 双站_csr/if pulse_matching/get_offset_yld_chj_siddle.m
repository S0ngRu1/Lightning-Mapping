yld_signal_start_loc = 469823486;
yld_signal_length = 1024;
yld_signal = read_signal('../20240822165932.6610CH1.dat',yld_signal_length,yld_signal_start_loc);
chj_signal_length = 1e8;
chj_signal_start_loc = 4.5e8;
chj_signal = read_signal('../2024 822 85933.651462CH1.dat',chj_signal_length,chj_signal_start_loc);
all_locs = [];
all_R_gccs = [];
all_t_gccs = [];
filtered_yld_signal = filter_bp(yld_signal, 20e6 ,80e6 ,5);
processed_yld_signal = real(windowsignal(detrend(filtered_yld_signal)));
subsignal_length = yld_signal_length/4;
subsignal_start = 1:subsignal_length:chj_signal_length;
h = waitbar(0, 'Processing...'); 
for subi = 1:numel(subsignal_start)
    waitbar(subi/numel(subsignal_start), h, sprintf('Processing %.2f%%', subi/numel(subsignal_start)*100));
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

[top10_r_gccs, top10_r_gcc_idxs] = maxk(all_R_gccs, 10);
start_read_loc_chj_top10 = zeros(1, 10);
for i = 1:10
    start_read_loc_chj_top10(i) = chj_signal_start_loc + subsignal_start(top10_r_gcc_idxs(i)) + floor(all_t_gccs(top10_r_gcc_idxs(i)));
end
close(h); 
offset_result = start_read_loc_chj_top10 - yld_signal_start_loc;
disp(offset_result)