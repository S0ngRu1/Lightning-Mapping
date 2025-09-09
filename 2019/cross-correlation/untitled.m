% 生成一个正弦信号
Fs = 1000; % 采样率
t = 0:1/Fs:1;
f = 50; % 信号频率
x = sin(2*pi*f*t);

% 进行傅里叶变换得到频域信号
X = fft(x);

% 在频域内应用窗函数
window = hamming(length(X)); % 使用汉明窗函数
X_windowed = X .* window';

% 进行逆傅里叶变换得到时域信号
x_windowed = ifft(X_windowed);

% 绘制原始信号和加窗处理后的信号对比
subplot(2,1,1);
plot(t, x);
title('原始信号');
xlabel('时间 (s)');
ylabel('幅度');
subplot(2,1,2);
plot(t, abs(x_windowed));
title('频域加窗处理后的信号');
xlabel('时间 (s)');
ylabel('幅度');
