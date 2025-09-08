
noise_std = 3.545;
signal_length = 200;
r_loction = 3.8e8+3500+9000+500;
ch1 = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction);
filtered_signal1 = filter_bp(ch1,30e6,80e6,5);
result = find_pulses_advanced(filtered_signal1,noise_std,200e6,3,5);

%% ========================================================================
%  可视化 find_pulses_advanced 算法的结果
% =========================================================================

fprintf('开始绘制算法结果...\n');

% --- 1. 准备绘图所需的数据 ---
waveform = filtered_signal1;
pulse_catalog = result;
sampling_rate_hz = 200e6;
ts_ns = 1 / sampling_rate_hz * 1e9; % 每个采样点的时间 (ns)

% 重新计算包络和阈值，用于绘图
envelope = abs(hilbert(waveform));
detection_threshold = noise_std * 3; % 与您调用函数时使用的因子保持一致
% 创建原始坐标点（采样点索引，从1开始）
sample_indices = 1:length(waveform);  % 原始采样点索引（1-based）

% --- 2. 创建图形并绘制基础信号 ---
figure;
hold on;

% 绘制滤波后的原始信号（使用采样点索引作为横轴）
plot(sample_indices, waveform, 'Color', [0.3, 0.7, 1.0], 'DisplayName', '预处理后信号');

% 绘制希尔伯特包络（原始坐标）
plot(sample_indices, envelope, 'r-', 'LineWidth', 1, 'DisplayName', '希尔伯特包络');

% 绘制检测阈值线（原始坐标）
plot(sample_indices, ones(size(sample_indices)) * detection_threshold, ...
    '--', 'Color', [1, 0.5, 0], 'LineWidth', 1, 'DisplayName', '检测阈值');

% --- 3. 在图上标记算法找到的脉冲 ---
if ~isempty(pulse_catalog)
    fprintf('正在标记 %d 个已识别的脉冲...\n', numel(pulse_catalog));
    for i = 1:numel(pulse_catalog)
        p = pulse_catalog(i);
        
        % --- 使用半透明灰色区域标记脉冲的起点和终点（直接用原始索引）---
        start_idx = p.start_idx;  % 原始起点索引
        end_idx = p.end_idx;      % 原始终点索引
        
        % 获取当前y轴范围，用于绘制区域
        y_lim = get(gca, 'YLim'); 
        
        % 用原始索引绘制脉冲区域（x轴为采样点索引）
        patch([start_idx, end_idx, end_idx, start_idx], ...
              [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], ...
              [0.5, 0.5, 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
              'HandleVisibility', 'off');
        
        % --- 用醒目的星号标记亚采样点精度的时刻（转换为原始坐标）---
        % 精确时刻的采样点坐标（可能为非整数，如3.2个采样点）
        precise_loc_samples = p.precise_time_ns / ts_ns + 1;  % 从时间转换回采样点索引
        precise_amplitude = interp1(sample_indices, envelope, precise_loc_samples);  % 用原始索引插值
        
        plot(precise_loc_samples, precise_amplitude, ...
             'p', 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'k', ...
             'MarkerSize', 5, 'HandleVisibility', 'off');
    end
    
    % 图例保持一致
    h_boundary = patch(NaN, NaN, [0.5, 0.5, 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h_timing = plot(NaN, NaN, 'p', 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'k', 'MarkerSize', 5);
    legend('滤波后信号', '希尔伯特包络', '检测阈值', '脉冲边界', '精确正时');
    
else
    fprintf('结果目录为空，无需标记。\n');
    legend('滤波后信号', '希尔伯特包络', '检测阈值');
end

% --- 4. 美化图形（更新坐标轴标签为原始坐标）---
hold off;
xlabel('采样点');  % 横轴改为采样点索引
ylabel('幅值');
set(gca, 'FontSize', 10);
xlim([min(sample_indices), max(sample_indices)]);  % 用原始索引限制x轴范围

fprintf('绘图完成！\n');