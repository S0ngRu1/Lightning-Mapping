fs = 200e6;
signal_length = 1024;
r_loction_yld = 4.698e8+23300;

% 读取信号
signal = read_signal('..\\20240822165932.6610CH1.dat', signal_length, r_loction_yld);
signal = filter_bp(signal,30e6,80e6,5);
% 构建滤波器组（默认 'amor'）
fb =cwtfilterbank( ...
    'SignalLength', length(signal), ...
    'SamplingFrequency', fs, ...
    'Wavelet', 'morse', ...
    'FrequencyLimits', [25e6 85e6], ...   % 强制覆盖更高频率
    'VoicesPerOctave', 48);              % 提高频率分辨率

% 进行小波变换
[cfs, f] = wt(fb, signal);

% 可视化
figure
scalogram_data = abs(cfs);
imagesc((1:length(signal))/fs*1e6, f/1e6, scalogram_data);  % us 和 MHz
axis xy;
xlabel('Time (\mus)');
ylabel('Frequency (MHz)');
title('Scalogram (CWT)');
colorbar;
