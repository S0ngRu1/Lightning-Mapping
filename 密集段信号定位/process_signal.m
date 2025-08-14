fs = 200e6; % 您的采样率
fp_start = 60e6; % 通带起始
fp_end = 75e6;   % 通带结束
fs_stop1 = 55e6; % 第一个阻带，要比50MHz更靠近通带，以获得陡峭滚降
fs_stop2 = 80e6; % 第二个阻带
signal_length = 10240;
r_loction_yld = 3.8e8+3500+6000;
ch1_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction_yld);
% 设计一个8阶巴特沃斯带通滤波器
order = 8;
[b, a] = butter(order, [fp_start, fp_end] / (fs/2), 'bandpass');


% 对密集段信号进行零相位滤波
processed_ch1_yld = filtfilt(b, a, ch1_yld);
plot(processed_ch1_yld);
plot_signal_spectrum(processed_ch1_yld); 


% -------------------------------------------------------------------------
% 用于计算能量包络的滑动窗口长度
energy_window_len = 50; 

min_peak_height = 10; % 峰值最小高度，
min_peak_distance = 128; 

% 以峰值为中心，进行处理的信号片段的总长度
processing_window_len = 256; 

fprintf('正在计算能量包络并寻找高质量事件...\n');
% 我们以CH1作为基准来寻找事件位置
% **修正了这里的变量名，使用您滤波后的信号**
energy_envelope = movmean(processed_ch1_yld.^2, energy_window_len);

% 寻找能量峰值的位置
[pks, locs] = findpeaks(energy_envelope, ...
                        'MinPeakHeight', min_peak_height, ...
                        'MinPeakDistance', min_peak_distance);

fprintf('找到了 %d 个高质量事件。\n', length(locs));

% --- 可视化检查，非常重要！---
% 绘制信号、能量包络和找到的峰值，检查参数设置是否合理
figure('Name', '高质量事件定位结果检查');
t_axis = (0:length(processed_ch1_yld)-1) / fs;
plot(t_axis, processed_ch1_yld, 'b-', 'DisplayName', '滤波后信号 (CH1)');
hold on;
plot(t_axis, energy_envelope, 'r-', 'LineWidth', 2, 'DisplayName', '能量包络');
plot(t_axis(locs), pks, 'kv', 'MarkerFaceColor', 'y', 'MarkerSize', 8, 'DisplayName', '找到的事件峰值');
legend('show');
title('高质量事件定位结果检查');
xlabel('时间 (s)');
ylabel('幅值 / 能量');
grid on;
zoom xon;


%% --- 进行精确加窗 ---
% -------------------------------------------------------------------------
fprintf('开始遍历事件，进行精确加窗和后续处理...\n');

if isempty(locs)
    fprintf('**警告：未找到任何满足条件的事件峰值，请降低 min_peak_height 参数！**\n');
    return;
end

% 创建一个汉宁窗
win = hann(processing_window_len);
win = win(:); % 确保是列向量

% 
figure('Name', '精确加窗效果展示');

for i = 1:length(locs)
    center_loc = locs(i);
    
    % 计算截取片段的起始和结束索引
    start_idx = center_loc - floor(processing_window_len / 2);
    end_idx = center_loc + ceil(processing_window_len / 2) - 1;
    
    % 边界检查，防止索引越界
    if start_idx < 1 || end_idx > length(processed_ch1_yld)
        fprintf('事件 %d (位置 %d) 过于靠近信号边界，已跳过。\n', i, center_loc);
        continue;
    end
    
    % 从三个通道截取信号片段 (确保是列向量)
    segment_ch1 = processed_ch1_yld(start_idx:end_idx);
    % 对截取的片段应用窗函数
    windowed_segment_ch1 = segment_ch1 .* win;

    if i == 1 
        subplot(2,1,1);
        plot(segment_ch1);
        title(sprintf('找到的第 %d 个事件片段 (未加窗)', i));
        grid on;
        
        subplot(2,1,2);
        plot(windowed_segment_ch1);
        title(sprintf('找到的第 %d 个事件片段 (已加窗)', i));
        grid on;
        xlabel('采样点');
        
    end
end
plot_signal_spectrum(windowed_segment_ch1); 
fprintf('所有高质量事件处理完毕。\n');


function signal = read_signal(signal_path, r_length,r_loction)
    fid  = fopen(signal_path,'r');%读取数据的位置

    %使用fseek函数将文件指针移动到指定位置，以便读取数据。
    %这里指定移动位置为r_location，表示移动到指定位置开始读取数据。
    fseek(fid,r_loction*2,'bof');
    %使用fread函数从文件中读取数据，读取的数据长度为r_length，数据以int16格式读取。
    %将读取到的数据分别保存到变量ch_1、ch_2和ch_3中。
    signal = fread(fid,r_length,'int16');
    %关闭所有文件
    fclose(fid);
end


function plot_signal_spectrum(signal)
% Plotting the signal spectrum时域信号的频谱图
fs = 200;
fft_signal = fft(signal);
n = length(fft_signal);
x = (0:n/2-1) * (fs/n);
figure
plot(x, 2.0 / n * abs(fft_signal(1:n/2)))

xlabel('Frequency (MHz)')
ylabel('Amplitude')
grid on
end



