fs = 200e6;
signal_length = 60;
r_loction_yld = 4.698e8 + 23300 + 720;

% 定义三个通道的文件名
channels = {'CH1', 'CH2', 'CH3'};
numChannels = length(channels);

% 为每个通道创建一个图形窗口，包含时频图和时域图
for i = 1:numChannels
    % 读取当前通道的信号
    filename = sprintf('..\\20240822165932.6610%s.dat', channels{i});
    signal = read_signal(filename, signal_length, r_loction_yld);
    signal = filter_bp(signal, 30e6, 80e6, 5);
    
    % 构建滤波器组（默认 'amor'）
    fb = cwtfilterbank( ...
        'SignalLength', length(signal), ...
        'SamplingFrequency', fs, ...
        'Wavelet', 'morse', ...
        'FrequencyLimits', [25e6 85e6], ...   % 强制覆盖更高频率
        'VoicesPerOctave', 48);              % 提高频率分辨率

    % 进行小波变换
    [cfs, f] = wt(fb, signal);
    
    % 创建一个新的图形窗口，包含两个子图
    figure('Name', sprintf('通道 %s 分析结果', channels{i}));
    
    % 绘制时频图
    subplot(2, 1, 1);
    scalogram_data = abs(cfs);
    imagesc((1:length(signal))/fs*1e6, f/1e6, scalogram_data);  % 单位转换为us和MHz
    axis xy;
    xlabel('时间 (\mus)');
    ylabel('频率 (MHz)');
    title(sprintf('通道 %s 时频图 (CWT)', channels{i}));
    colorbar;
    
    % 绘制时域波形
    subplot(2, 1, 2);
    plot((1:length(signal))/fs*1e6, signal);  % 时间轴单位为us
    xlabel('时间 (\mus)');
    ylabel('信号幅度');
    title(sprintf('通道 %s 时域波形', channels{i}));
    grid on;
end
