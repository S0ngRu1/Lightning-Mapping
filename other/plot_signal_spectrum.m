function plot_signal_spectrum(signal)
% Plotting the signal spectrum时域信号的频谱图
fs = 400;
fft_signal = fft(signal);
n = length(fft_signal);
x = (0:n/4-1) * (fs/n);
figure
plot(x, 4.0 / n * abs(fft_signal(1:n/4)))

xlabel('Frequency (MHz)')
ylabel('Amplitude')
grid on
end