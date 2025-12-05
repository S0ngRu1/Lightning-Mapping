%% 在一幅图中绘制通道1和通道4的图

signal_length = 5e5;
r_loction_yld = 369650000;
ch4_yld = read_signal('..\\20240822165932.6610CH4.dat', signal_length, r_loction_yld);
ch1_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction_yld);
bp_filtered_ch1 = filter_bp(ch1_yld,30e6,80e6,5);
x_indices = 0 :  signal_length - 1;
time_ms = x_indices * (5 / 1e3);
baseline = movmedian(ch4_yld, 1024);
E_fast = ch4_yld - baseline;
bp_filtered_E_fast = filter_bp(E_fast,10e6,90e6,5);
figure
title(sprintf('%d + %d 电场', r_loction_yld, signal_length));
% plot(bp_filtered_ch1+300)
% hold on 
% plot(bp_filtered_E_fast)
plot(time_ms,bp_filtered_ch1+300)
hold on 
plot(time_ms,bp_filtered_E_fast)
xlabel('时间 (us)');
