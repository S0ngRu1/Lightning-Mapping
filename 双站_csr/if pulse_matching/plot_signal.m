signal_length = 3e7; 
start_read_loc_yld = 3.7e8;
start_read_loc_chj = start_read_loc_yld+34151156;
yld_signal = read_signal('../20240822165932.6610CH1.dat',signal_length,start_read_loc_yld);
chj_signal = read_signal('../2024 822 85933.651462CH1.dat',signal_length,start_read_loc_chj);
x = downsample(yld_signal,1);
y = downsample(chj_signal,1);
% plot(x);
subplot(2,1,1);plot(x);title('yld');xlabel('采样点');ylabel('振幅');
subplot(2,1,2);plot(y);title('chj');xlabel('采样点');ylabel('振幅');

filtered_yld_signal1 = filter_bp(yld_signal,30e6,80e6,5);
filtered_chj_signal1 = filter_bp(chj_signal,30e6,80e6,5);
subplot(2,1,1);plot(filtered_yld_signal1);title('filtered_yld');xlabel('采样点');ylabel('振幅');
subplot(2,1,2);plot(filtered_chj_signal1);title('filtered_chj');xlabel('采样点');ylabel('振幅');