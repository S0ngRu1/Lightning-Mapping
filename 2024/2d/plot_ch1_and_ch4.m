%% 在一幅图中绘制通道1和通道4的图

signal_length = 0.5e8;
r_loction_yld = 3.65e8+1.2e8;
ch4_yld = read_signal('..\\20240822165932.6610CH4.dat', signal_length, r_loction_yld);
ch1_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction_yld);
bp_filtered_ch1 = filter_bp(ch1_yld,30e6,80e6,5);
x_indices = 0 :  signal_length - 1;
time_ms = x_indices * (5 / 1e3);
baseline = movmedian(ch4_yld, 1024);
E_fast = ch4_yld - baseline;
figure
% plot(bp_filtered_ch1+300)
% hold on 
% plot(E_fast)
plot(time_ms,bp_filtered_ch1+300)
hold on 
plot(time_ms,E_fast)
title(sprintf('%d + %d 电场', r_loction_yld, signal_length));
xlabel('时间 (us)');
