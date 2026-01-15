clear; clc; close all;

%% ==================== 1. 参数设置 ====================
% 文件路径 (根据您的描述修改)
raw_file = '..\20240822165932.6610CH1.dat'; 

% 采样率
Fs = 200e6; 

% --- 【关键】：请根据之前的定位结果，挑选一段典型的信号位置 ---
% 单位：采样点 (Sample Index)
% 建议选择一段 200us ~ 500us 长的数据 (约 40000 ~ 100000 点)
Start_Loc_Neg = 367000000+2800000; % [示例] 负先导大概位置 (Pre-RS)
Start_Loc_Pos = 383000000+1400000; % [示例] 正先导大概位置 (Post-RS)

% 读取长度 (例如 200 us)
T_duration_us = 200; 
Read_Len = round(T_duration_us * 1e-6 * Fs);

%% ==================== 2. 数据读取与预处理 ====================
% 2.1 读取负先导
if ~isfile(raw_file), error('找不到文件: %s', raw_file); end
sig_neg_raw = read_signal(raw_file, Read_Len, Start_Loc_Neg);
sig_neg_filt = filter_bp(sig_neg_raw, 30e6, 80e6, Fs, 5); % 30-80MHz 滤波

% 2.2 读取正先导
sig_pos_raw = read_signal(raw_file, Read_Len, Start_Loc_Pos);
sig_pos_filt = filter_bp(sig_pos_raw, 30e6, 80e6, Fs, 5); 

% 2.3 构建时间轴 (微秒)
t_axis = (0:length(sig_neg_filt)-1) / Fs * 1e6;

%% ==================== 3. 频谱分析 (FFT) ====================
% 计算负先导频谱
[f_axis, amp_neg] = calc_spectrum(sig_neg_filt, Fs);

% 计算正先导频谱
[f_axis, amp_pos] = calc_spectrum(sig_pos_filt, Fs);

%% ==================== 4. JGR 风格绘图 ====================
figure('Units', 'centimeters', 'Position', [5, 5, 18, 14], 'Color', 'w');
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 字体设置
font_name = 'Arial';
font_size = 9;
label_size = 10;

% --- (a) 负先导时域波形 ---
nexttile;
plot(t_axis, sig_neg_filt, 'b', 'LineWidth', 0.5);
xlim([0, T_duration_us]);
% 为了美观，自动调整Y轴范围并居中
ymax = max(abs(sig_neg_filt)) * 1.1; ylim([-ymax, ymax]);
ylabel('Amplitude', 'FontName', font_name, 'FontSize', label_size);
title('(a) Negative Leader Waveform (Time Domain)', 'FontName', font_name, 'FontSize', label_size, 'FontWeight', 'bold');
apply_jgr_style(gca, font_name, font_size);

% --- (b) 正先导时域波形 ---
nexttile;
plot(t_axis, sig_pos_filt, 'r', 'LineWidth', 0.5);
xlim([0, T_duration_us]);
ymax = max(abs(sig_pos_filt)) * 1.1; ylim([-ymax, ymax]);
% ylabel('Amplitude (ADU)'); % 右侧省略Y轴标签
title('(b) Positive Leader Waveform (Time Domain)', 'FontName', font_name, 'FontSize', label_size, 'FontWeight', 'bold');
apply_jgr_style(gca, font_name, font_size);

% --- (c) 负先导频谱 ---
nexttile;
plot(f_axis, amp_neg, 'b', 'LineWidth', 1.0);
xlim([0, 100]); % 关注 0-100 MHz
xlabel('Frequency (MHz)', 'FontName', font_name, 'FontSize', label_size);
ylabel('Normalized Amplitude', 'FontName', font_name, 'FontSize', label_size);
title('(c) Negative Leader Spectrum', 'FontName', font_name, 'FontSize', label_size, 'FontWeight', 'bold');
apply_jgr_style(gca, font_name, font_size);
% 标注峰值区域
x_fill = [60 75 75 60]; y_fill = [0 0 max(amp_neg) max(amp_neg)];
% patch(x_fill, y_fill, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 可选：高亮 60-75MHz

% --- (d) 正先导频谱 ---
nexttile;
plot(f_axis, amp_pos, 'r', 'LineWidth', 1.0);
xlim([0, 100]);
xlabel('Frequency (MHz)', 'FontName', font_name, 'FontSize', label_size);
% ylabel('Normalized Amplitude'); 
title('(d) Positive Leader Spectrum', 'FontName', font_name, 'FontSize', label_size, 'FontWeight', 'bold');
apply_jgr_style(gca, font_name, font_size);

% 增加整体标题 (可选)
% sgtitle('Comparison of Radiation Signals: Negative vs. Positive Leader', 'FontName', font_name, 'FontSize', 12);

%% ==================== 5. 辅助函数 ====================

function apply_jgr_style(ax, fname, fsize)
    set(ax, 'FontName', fname, 'FontSize', fsize, 'LineWidth', 1.0, ...
        'Box', 'on', 'TickDir', 'in', 'GridAlpha', 0.3, 'GridLineStyle', ':');
    grid on;
end

function [f, P1] = calc_spectrum(signal, Fs)
    % 计算单边幅度谱
    L = length(signal);
    
    % 加窗 (汉明窗) 以减少频谱泄漏 - 这里的窗加在时域
    w = hamming(L);
    signal_win = signal .* w;
    
    Y = fft(signal_win);
    
    P2 = abs(Y/L);       % 双边谱
    P1 = P2(1:floor(L/2)+1);    % 单边谱
    P1(2:end-1) = 2*P1(2:end-1);
    
    % 归一化 (方便对比形态，看相对分布)
    P1 = P1 / max(P1);
    
    f = Fs*(0:(L/2))/L;  % 频率轴 (Hz)
    f = f / 1e6;         % 转换为 MHz
end

function filtered_signal = filter_bp(signal, f1, f2, Fs, order)
    % 巴特沃斯带通滤波器
    fn = Fs/2;
    Wn = [f1 f2]/fn;
    [b,a] = butter(order, Wn); 
    filtered_signal = filtfilt(b, a, double(signal)); % 转换为double防止溢出
end

function signal = read_signal(signal_path, r_length, r_location)
    fid = fopen(signal_path, 'r');
    if fid == -1, error('无法打开文件: %s', signal_path); end
    % fseek 移动文件指针
    status = fseek(fid, r_location * 2, 'bof'); % *2 因为 int16 是 2 字节
    if status == -1, fclose(fid); error('fseek 失败'); end
    % 读取数据
    signal = fread(fid, r_length, 'int16');
    fclose(fid);
    
    % 简单的去零漂 (减去均值)
    signal = signal - mean(signal);
end