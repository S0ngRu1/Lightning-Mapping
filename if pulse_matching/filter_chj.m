signal_length = 6000;
start_read_loc_chj = 4e8;
chj1 = read_signal('../2024 822 85933.651462CH1.dat',signal_length,start_read_loc_chj);
chj2 = read_signal('../2024 822 85933.651462CH2.dat',signal_length,start_read_loc_chj);
chj3 = read_signal('../2024 822 85933.651462CH3.dat',signal_length,start_read_loc_chj+165/5);
filtered_chj3_signal = filter_bp(chj3, 20e6, 80e6, 5);
figure;
plot(chj3)
% 对信号进行滤波
chj3 = filter(ones(1, 5)/5, 1, chj3);

figure;
plot(chj3)

figure;
plot(filtered_chj3_signal)

