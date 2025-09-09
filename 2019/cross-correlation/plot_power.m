function plot_power(signal)
Y =fft(detrend(signal));       %快速FFT变换
N = length(Y);    %FFT变换后数据长度
power = abs(Y(1:N/2)).^2;  %求功率谱
nyquist = 200e6/2;
freq = (1:N/2)/(N/2)*nyquist; %求频率
plot(freq,power); grid on     %绘制功率谱图
xlabel('频率')
ylabel('功率')
title('功率谱图')
end