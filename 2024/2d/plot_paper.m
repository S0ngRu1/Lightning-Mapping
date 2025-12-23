%% ========================================================================
%  创新点1专项论证：自适应窗口捕获机制
% =========================================================================
clear; clc; close all;
% --- 1. 参数准备与长信号读取 ---
fs = 200e6;
ts = 1/fs;
block_len = 2e6;
r_location = 4.02e8; % 建议选一个脉冲密集的起始点
% 请确保路径正确
ch1 = read_signal('..\\20240822165932.6610CH1.dat', block_len, r_location);
filtered_sig = filter_bp(detrend(ch1), 30e6, 80e6, 5);
envelope = abs(hilbert(filtered_sig));
% --- 2. 模拟/计算密度分布与自适应长度 ---
% 定义滑动窗口来统计脉冲密度
density_win = 1e5;
num_sections = floor(block_len / density_win);
pulse_counts = zeros(1, num_sections);
window_lengths = zeros(1, num_sections);
for i = 1:num_sections
    idx = (i-1)*density_win + (1:density_win);
    p_indices = find_pulses_advanced(filtered_sig(idx), 3.545, fs, 6, 5);
    cnt = numel(p_indices);
    pulse_counts(i) = cnt;
    window_lengths(i) = simulate_adaptive_logic(cnt);
end

% 时间轴生成
time_ms = (0:block_len-1) * ts * 1e3;
time_sections = (0:num_sections-1) * (density_win * ts * 1e3) + (density_win/2*ts*1e3);


ch1_A = read_signal('..\\20240822165932.6610CH1.dat',512,3.8e8+3500+9000+1024+1024);
filtered_signal_A = filter_bp(ch1_A,30e6,80e6,5);
pulse_catalog_micro = find_pulses_advanced(filtered_signal_A,3.2,200e6,4,10);
% --- 1. 准备绘图所需的数据 ---
waveform_A = filtered_signal_A;
sampling_rate_hz = 200e6;
ts_ns = 1 / sampling_rate_hz * 1e9;
envelope_A = abs(hilbert(waveform_A));
detection_threshold_A = 3 * 4;
sample_indices_A = 1:length(waveform_A);

% --- 3. 开始绘图 ---
figure('Color', 'w', 'Units', 'pixels', 'Position', [100, 50, 900, 1000]);
t = tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'tight');

% -------------------------------------------------------------------------
% 子图 1: 微观脉冲识别原理 
% -------------------------------------------------------------------------
ax1 = nexttile; % 激活第1块
hold(ax1, 'on'); box(ax1, 'on'); grid(ax1, 'on');
set(ax1, 'LineWidth', 1.2);

% [层级1] 噪声基底
fill([sample_indices_A, fliplr(sample_indices_A)], ...
    [ones(size(sample_indices_A))*-detection_threshold_A, ones(size(sample_indices_A))*detection_threshold_A], ...
    [0.95, 0.95, 0.95], 'EdgeColor', 'none', 'DisplayName', 'Noise Floor');
% [层级2] 阈值线
yline(detection_threshold_A, '--', 'Color', [0.85, 0.33, 0.1], 'LineWidth', 1.5, ...
    'DisplayName', 'Threshold (3\sigma)');
% [层级3] 波形与包络
p_wave = plot(sample_indices_A, waveform_A, '-', 'Color', [0.2, 0.6, 0.8, 0.6], 'LineWidth', 1, 'DisplayName', 'Waveform');
p_env = plot(sample_indices_A, envelope_A, '-', 'Color', [0.8, 0.1, 0.1], 'LineWidth', 1.5, 'DisplayName', 'Envelope');

% [层级4] 脉冲标记
y_limits = max(abs(waveform_A)) * 1.3;
ylim(ax1, [-y_limits, y_limits]);

if ~isempty(pulse_catalog_micro)
    for i = 1:numel(pulse_catalog_micro)
        p = pulse_catalog_micro(i);
        p_start_t = sample_indices_A(p.start_idx);
        p_end_t = sample_indices_A(p.end_idx);
        p_peak_t = sample_indices_A(p.peak_loc);
        
        % A. 绿色区域
        fill([p_start_t, p_end_t, p_end_t, p_start_t], ...
            [-y_limits, -y_limits, y_limits, y_limits], ...
            [0.2, 0.8, 0.2], 'FaceAlpha', 0.15, 'EdgeColor', 'none'); 
        % B. 五角星
        plot(p_peak_t, envelope_A(p.peak_loc), 'p', ...
            'MarkerFaceColor', [1, 0.8, 0], 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
        % C. 文字
        text(p_peak_t, envelope_A(p.peak_loc) + y_limits*0.1, ...
            ['P', num2str(i)], 'HorizontalAlignment', 'center', ...
            'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.8, 0.1, 0.1]);
    end
end

% 状态框 (固定在左上角)
str = {'\bf 算法状态分析', ...
    ['检测脉冲数: \color{blue}' num2str(numel(pulse_catalog_micro))], ...
    '当前密度: \color{red} 高 (Burst)', ...
    '窗口策略: \bf 缩短窗口'};

% 坐标解释：
% 0.98: 靠右侧 (1.0 是最右边)
% 0.03: 靠底部 (0.0 是最底部)
% HorizontalAlignment: 'right' -> 文本框向左生长
% VerticalAlignment: 'bottom' -> 文本框向上生长
text(ax1, 0.98, 0.03, str, 'Units', 'normalized', ...
    'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontSize', 9, ...
    'HorizontalAlignment', 'right', ... 
    'VerticalAlignment', 'bottom');

title('图A：脉冲发现算法效果', 'FontWeight', 'bold');
ylabel('幅值');
xlabel('采样点');
xlim([min(sample_indices_A), max(sample_indices_A)]);
legend([p_wave, p_env], '原始波形', '希尔伯特包络', 'Location', 'northeast');

% -------------------------------------------------------------------------
% 子图 2: 宏观脉冲密度 (Macro Density)
% -------------------------------------------------------------------------
ax2 = nexttile; % 激活第2块
hold(ax2, 'on'); grid(ax2, 'on');

% --- 左轴：绘制原始波形 + 包络 ---
yyaxis left
% 1. [新增] 绘制原始波形 (淡蓝色，细线，作为背景)
% 注意：变量名需与前面定义保持一致 (filtered_macro)
plot(time_ms, filtered_sig, '-', 'Color', [0.6, 0.8, 1, 0.8], ...
    'LineWidth', 0.5, 'DisplayName', '原始波形');

% 2. 绘制信号包络 (深灰色，叠加在波形之上)
plot(time_ms, envelope, '-', 'Color', [0.4, 0.4, 0.4], ...
    'LineWidth', 0.5, 'DisplayName', '信号包络');

ylabel('幅值');
set(gca, 'YColor', 'k'); % 左轴颜色设为黑色

% --- 右轴：绘制脉冲密度 ---
yyaxis right
plot(time_sections, pulse_counts, '-o', 'Color', [0.85, 0.33, 0.1], ...
    'LineWidth', 1.5, 'DisplayName', '脉冲密度');
ylabel('脉冲计数 / 窗口');
set(gca, 'YColor', [0.85, 0.33, 0.1]); % 右轴颜色设为橙色

% --- 修饰与标注 ---
xlabel('时间 (ms)');
xlim([min(time_ms), max(time_ms)]);
title('图B：宏观脉冲密度分布', 'FontWeight', 'bold');

% 标注高密度区
[max_cnt, max_idx] = max(pulse_counts);
if max_idx <= length(time_sections)
    % 这里的 y 坐标取包络的最大值，避免文字与波形重叠
    text(time_sections(max_idx), max(envelope)*0.7, '\leftarrow 脉冲密集区', ...
        'FontWeight', 'bold', 'Color', 'r', 'FontSize', 10);
end

% 更新图例 (包含原始波形)
legend('Location', 'northeast');

% -------------------------------------------------------------------------
% 子图 3: 对比展示 (Comparison)
% -------------------------------------------------------------------------
ax3 = nexttile; % 激活第3块
hold(ax3, 'on'); grid(ax3, 'on');
zoom_start = (max_idx-1)*density_win + 4000;

zoom_len = 1200;

zoom_idx = zoom_start : (zoom_start + zoom_len);

% 防止索引越界

if zoom_idx(end) > length(time_ms)

zoom_idx = length(time_ms)-zoom_len : length(time_ms);

end

time_zoom = time_ms(zoom_idx);

sig_zoom = filtered_sig(zoom_idx);

plot(time_zoom, sig_zoom, 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1);

% 1. 模拟固定窗口 (红色虚线 - 强行切断)

fixed_len = 500;

% 确保绘图框不越界

if fixed_len < length(time_zoom)

fixed_box_x = [time_zoom(1), time_zoom(fixed_len), time_zoom(fixed_len), time_zoom(1)];

patch(fixed_box_x, [-max(sig_zoom) -max(sig_zoom) max(sig_zoom) max(sig_zoom)], ...
'r', 'FaceAlpha', 0.1, 'EdgeColor', 'r', 'LineStyle', '--', 'LineWidth', 1.2);

text(time_zoom(fixed_len), -max(sig_zoom)*0.8, ' \times 信号截断 (固定窗口)', ...
'Color', 'r', 'FontWeight', 'bold', 'FontSize', 9);

end

% 2. 模拟自适应窗口 (蓝色实线 - 完美包裹)

adaptive_len = 1000;

if adaptive_len < length(time_zoom)

adapt_box_x = [time_zoom(1), time_zoom(adaptive_len), time_zoom(adaptive_len), time_zoom(1)];

patch(adapt_box_x, [-max(sig_zoom)*1.1 -max(sig_zoom)*1.1 max(sig_zoom)*1.1 max(sig_zoom)*1.1], ...
'b', 'FaceAlpha', 0.05, 'EdgeColor', 'b', 'LineWidth', 1.5);

text(time_zoom(adaptive_len), max(sig_zoom)*0.8, ' \surd 完整保留 (自适应)', ...
'Color', 'b', 'FontWeight', 'bold', 'FontSize', 9);

end

xlabel('时间 (ms)');

ylabel('归一化幅值');

title('图C：自适应窗口与固定窗口对比 (微观细节)', 'FontWeight', 'bold');

xlim([min(time_zoom) max(time_zoom)]);
% -------------------------------------------------------------------------
% 子图 4: 逻辑阶梯图 (Logic Staircase)
% -------------------------------------------------------------------------
ax4 = nexttile; % 激活第4块
hold(ax4, 'on'); grid(ax4, 'on');
% --- 使用 bar 绘制连续柱状图 ---
% 第三个参数 1 表示 BarWidth=1，让柱子紧密相连，形成“从这一点到下一点”的视觉效果
b = bar(time_sections, window_lengths, 1);

% --- 美化柱子外观 ---
b.FaceColor = [0.3, 0.7, 0.4]; % 设置为清新的绿色
b.EdgeColor = 'none';          % [关键] 去掉柱子边框。如果数据很密，保留边框会让图变成一团黑
b.FaceAlpha = 0.8;             % 设置一点透明度，看起来更现代

% --- 坐标轴设置 ---
ylabel('窗口长度');
xlabel('时间 (ms)');
title('图D：自适应窗口长度调节机制', 'FontWeight', 'bold');

% 固定 Y 轴刻度，只显示你的这几档窗口
ylim([0, 4500]);
yticks([512, 1024, 2048, 4096]);

% 确保 X 轴范围与其他子图对齐
% 注意：这里使用 time_ms_macro 的范围，以保证和上面几个图对齐
xlim([min(time_ms), max(time_ms)]);
% --- 辅助函数 ---
function w_len = simulate_adaptive_logic(count)
    if count < 5, w_len = 4096;
    elseif count < 15, w_len =  2048;
    elseif count < 25, w_len = 1024;
    else, w_len = 512;
    end
end