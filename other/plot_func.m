signal_length = 3e5;
r_loction = 7e8;
ch1 = read_signal_tdms('20250820151326_1505CH1.tdms',signal_length,r_loction);
figure
plot(ch1)
plot_signal_spectrum(ch1)