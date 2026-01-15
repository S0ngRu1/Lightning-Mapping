clear; clc; close all;

%% === 1. 参数设置 ===
filename = '..\2024\2d\results\result_yld_3.65e8_5.6e8_window_ADAPTIVE_1e4_factor3_with_error.txt';
fs = 200e6;
Start_loc_Base = 469200000;
Error_Threshold = 0.5;
Rcorr_Default = 0.3;

% --- 【关键：设置局部放大区域 (第四个子图)】 ---
% 请根据您的数据分布，手动调整下面这个范围，
% 选一个通道分叉或者拐弯的地方，范围大约跨度 20-30度
Zoom_Az_Lim = [160, 190];  % 方位角放大范围
Zoom_El_Lim = [30, 50];    % 仰角放大范围

%% === 2. 读取与筛选 (保持不变) ===
if ~isfile(filename), error('文件不存在'); end
opts = detectImportOptions(filename); opts.VariableNamingRule = 'preserve';
T = readtable(filename, opts);

valid_geo = T.Start_loc > Start_loc_Base & T.Start_loc < (Start_loc_Base + 1800000) & ...
            T.Elevation < 85 & T.Azimuth < 250 & abs(T.t123) < 10.0;
bad_region = false(height(T), 1);
if ismember('Azimuth', T.Properties.VariableNames)
    bad_region = (T.Azimuth > 160) & (T.Azimuth < 250) & (T.Elevation > 0) & (T.Elevation < 35);
end

mask_rcorr = false(height(T), 1);
if ismember('Win_len', T.Properties.VariableNames)
    idx_512 = (T.Win_len == 512);   mask_rcorr(idx_512) = T.Rcorr(idx_512) > 0.3;
    idx_1024 = (T.Win_len == 1024); mask_rcorr(idx_1024) = T.Rcorr(idx_1024) > 0.1;
    idx_2048 = (T.Win_len == 2048); mask_rcorr(idx_2048) = T.Rcorr(idx_2048) > 0.1;
    idx_4096 = (T.Win_len == 4096); mask_rcorr(idx_4096) = T.Rcorr(idx_4096) > 0.1;
    other = ~ismember(T.Win_len, [512, 1024, 2048, 4096]);
    mask_rcorr(other) = T.Rcorr(other) > Rcorr_Default;
else
    mask_rcorr = T.Rcorr > Rcorr_Default;
end

mask_error = (T.Err_Az < Error_Threshold) & (T.Err_El < Error_Threshold);
final_idx = valid_geo & (~bad_region) & mask_rcorr & mask_error;
data = T(final_idx, :);
data.Time_us = (data.Start_loc - Start_loc_Base) / fs * 1e6;

%% === 3. 绘图 (2x2 布局) ===
% 设置宽大于长的尺寸 (单位: 厘米)，适合双栏排版
fig_width = 18;  % 宽
fig_height = 14; % 高
f = figure('Units', 'centimeters', 'Position', [5, 5, fig_width, fig_height], 'Color', 'w');

% JGR 字体风格
font_name = 'Arial';
font_size = 10;
label_size = 11;

% 2行2列 紧凑布局
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
colors = data.Time_us;
time_range = [min(data.Time_us), max(data.Time_us)];

% --- (a) 方位角误差 (左上) ---
ax1 = nexttile;
scatter(data.Time_us, data.Err_Az, 8, colors, 'filled', 'MarkerFaceAlpha', 0.6);
hold on; grid on;
ylabel('Azimuth Error (°)', 'FontName', font_name, 'FontSize', label_size);
xlim(time_range); ylim([0, Error_Threshold]);
apply_jgr_style(ax1, font_name, font_size, '(a)');
xticklabels([]); % 移除X轴标签

% --- (b) 仰角误差 (左下) ---
ax2 = nexttile;
scatter(data.Time_us, data.Err_El, 8, colors, 'filled', 'MarkerFaceAlpha', 0.6);
hold on; grid on;
xlabel('Time (\mus)', 'FontName', font_name, 'FontSize', label_size);
ylabel('Elevation Error (°)', 'FontName', font_name, 'FontSize', label_size);
xlim(time_range); ylim([0, Error_Threshold]);
apply_jgr_style(ax2, font_name, font_size, '(b)');

% --- (c) 二维全景图 (右上) ---
ax3 = nexttile;
hold on; axis equal; grid on;
% 绘制灰色误差棒 (细一点)
errorbar(data.Azimuth, data.Elevation, data.Err_El, data.Err_El, data.Err_Az, data.Err_Az, ...
    '.', 'Color', [0.8 0.8 0.8], 'CapSize', 0, 'LineWidth', 0.1); 
scatter(data.Azimuth, data.Elevation, 4, colors, 'filled'); % 点小一点
ylabel('Elevation (°)', 'FontName', font_name, 'FontSize', label_size);
xticklabels([]); % 移除X轴标签，保持整洁
apply_jgr_style(ax3, font_name, font_size, '(c) Overview');

% 在全景图上画红框，标示放大区域
rectangle('Position', [Zoom_Az_Lim(1), Zoom_El_Lim(1), diff(Zoom_Az_Lim), diff(Zoom_El_Lim)], ...
          'EdgeColor', 'r', 'LineWidth', 1.2, 'LineStyle', '-');
% 自动调整全景图范围
margin = 2;
xlim([min(data.Azimuth)-margin, max(data.Azimuth)+margin]);
ylim([min(data.Elevation)-margin, max(data.Elevation)+margin]);

% --- (d) 局部放大图 (右下) ---
ax4 = nexttile;
hold on; axis equal; grid on;
% 1. 先画误差棒 (颜色深一点 [0.6 0.6 0.6]，显眼一点)
eb = errorbar(data.Azimuth, data.Elevation, data.Err_El, data.Err_El, data.Err_Az, data.Err_Az, ...
    '.', 'Color', [0.6 0.6 0.6], 'CapSize', 0, 'LineWidth', 0.8);
eb.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 2. 再画数据点 (点稍微大一点，突出位置)
scatter(data.Azimuth, data.Elevation, 15, colors, 'filled');

% 3. 设置强制放大范围
xlim(Zoom_Az_Lim);
ylim(Zoom_El_Lim);

xlabel('Azimuth (°)', 'FontName', font_name, 'FontSize', label_size);
ylabel('Elevation (°)', 'FontName', font_name, 'FontSize', label_size);
apply_jgr_style(ax4, font_name, font_size, '(d) Zoomed View');

% 设置边框颜色为红色，呼应上一张图的红框
set(ax4, 'XColor', 'r', 'YColor', 'r', 'LineWidth', 1.2); 

% --- 共享 Colorbar ---
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Time (\mus)';
cb.Label.FontName = font_name; cb.Label.FontSize = label_size;
colormap('jet');

% 链接左侧时间轴
linkaxes([ax1, ax2], 'x');

%% 辅助函数
function apply_jgr_style(ax, fname, fsize, label_text)
    set(ax, 'FontName', fname, 'FontSize', fsize, 'LineWidth', 1.0, ...
        'Box', 'on', 'TickDir', 'in', 'GridAlpha', 0.3, 'GridLineStyle', ':');
    % 标签放在左上角
    text(ax, 0.03, 0.93, label_text, 'Units', 'normalized', ...
        'FontName', fname, 'FontSize', fsize+1, 'FontWeight', 'bold');
end