signal_length = 1024; 
yld_signal = read_signal('../20240822165932.6610CH1.dat',signal_length,470799897);
chj_signal = read_signal('../2024 822 85933.651462CH1.dat',signal_length,504998496);
x = downsample(yld_signal,1);
y = downsample(chj_signal,1);
% plot(x);
subplot(2,1,1);plot(x);title('yld');xlabel('采样点');ylabel('振幅');
subplot(2,1,2);plot(y);title('chj');xlabel('采样点');ylabel('振幅');
