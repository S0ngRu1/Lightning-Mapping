% --- 0. 清理环境 ---
clear; clc; close all;


% --- 2. 准备数据 ---
sampling_points = 50;
original_x = 1:sampling_points;
peak_center = 25.3; 
peak_width = 4;
original_signal = exp(-((original_x - peak_center).^2) / (2 * peak_width^2));
original_signal = original_signal + randn(1, sampling_points) * 0.03;

% 设置上采样倍率
upsampling_factor = 20;

% 调用函数进行上采样
upsampled_data = upsampling(original_signal, upsampling_factor);
x_up = upsampled_data(1, :);
signal_up = upsampled_data(2, :);

% --- 3. 绘制最终对比图 ---
figure('Position', [100, 100, 1000, 600]);
hold on;

% 步骤 A: 先绘制平滑的插值曲线
plot(x_up, signal_up, '-', ...
    'Color', [0.8500, 0.3250, 0.0980], ... % 橙色
    'LineWidth', 2.5);

% 步骤 B: 再绘制原始的、带连线的离散采样点
% << 主要修改点：将 'o' 改为 'o-'，表示同时绘制圆圈和连线
plot(original_x, original_signal, 'o-', ... 
    'Color', [0, 0.4470, 0.7410], ... % 蓝色
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'w', ...
    'LineWidth', 1.5); % 适当调整线宽以区分

hold off;

% --- 4. 美化和标注 ---
title('上采样前后信号对比 ', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('采样点索引', 'FontSize', 12);
ylabel('信号幅值', 'FontSize', 12);
grid on;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

% 定义要放大的区域
zoom_center_idx = round(peak_center);
zoom_width = 4;
xlim([zoom_center_idx - zoom_width, zoom_center_idx + zoom_width]);

% << 更新图例说明
legend('上采样后的平滑曲线', ...
       '原始信号 ', ...
       'Location', 'northwest', 'FontSize', 11);

% --- 1. 使用您原来的 upsampling 函数 ---
% 这个函数返回一个 2xN 的矩阵，第一行是x，第二行是y
function new_signal = upsampling(original_signal, upsampling_factor)
    original_x = (1:numel(original_signal))';
    original_y = original_signal;
    upsampled_length = length(original_x) * upsampling_factor;
    upsampled_x = linspace(1, length(original_x), upsampled_length);
    interpolated_signal = interp1(original_x, original_y, upsampled_x, 'spline');
    % 返回一个 2xN 的矩阵
    new_signal = [upsampled_x; interpolated_signal]; 
end
