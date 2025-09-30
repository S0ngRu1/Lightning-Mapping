signal_length = 1e4;
r_loction_yld = 4.698e8;
ch1_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction_yld);
bp_filtered_yld = filter_bp(ch1_yld,20e6,80e6,5);
[high_freq_signal, dominant_freqs] = extract_band_emd(bp_filtered_yld, 200e6, 10e6, 100e6, true);
