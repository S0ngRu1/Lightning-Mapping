%读取数据
signal_length = 2e3;
sig = read_signal('20190604164852.7960CH1.dat',signal_length);
fs = 200e6;
t = (0:signal_length-1) * 1e6/fs; % 将时间单位从秒改为微秒
figure('color','white')
plot(t, sig, 'k') % 绘制原始信号
xlabel('Time (\mu s)') % 设置 x 轴标签单位为微秒
ylabel('Amplitude')
title('Original Signal')
%EEMD分解
Nstd = 0.2; %Nstd为附加噪声标准差与Y标准差之比
NE = 20;   %NE为对信号的平均次数
imf = eemd(sig,Nstd,NE);

%绘制EEMD的各个IMF图和其对应的频谱图
i = size(imf, 2);
figure;
for j = 1:i
    % 画出IMF信号图
    subplot(i, 1, j);
    plot(imf(:,j));
    title(['IMF' num2str(j) ' 信号']);
    
end

figure;
for j = 1:i
% 计算并绘制IMF信号的频谱图
    subplot(i, 1, j);
    fs = 200e6;
    fft_signal = fft(imf(:,j));
    n = length(fft_signal);
    x = (0:n/2-1) * (fs/n);
    plot(x, 2.0 / n * abs(fft_signal(1:n/2)));
    title(['IMF' num2str(j) ' 频谱']);
end


%信号重构
indices = [3];  
filtered_imfs = imf(:,indices);
filtered_signal1 = sum(filtered_imfs, 2);
ori = sig;  %无噪声信号
fil = filtered_signal1;  %滤波后信号
figure('color','w')
subplot(211);plot(sig,'k');title('原始信号')
subplot(212);plot(fil,'k');title('滤波后信号')

%带通滤波
filtered_bandpass = bandpass(fil, [40e6 80e6],200e6);


plot_signal_spectrum(fil);

plot_signal_spectrum(sig);

plot_signal_spectrum(filtered_bandpass);

figure('color','w')
subplot(211);plot(fil,'k');title('原始信号')
subplot(212);plot(filtered_bandpass,'k');title('滤波后信号')

[r12,lags12] = xcorr(sig,filtered_bandpass,'normalized');

figure
plot(lags12,r12)
