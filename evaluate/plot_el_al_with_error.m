clear; clc; close all;

%% === 1. 全局参数与文件设置 (必须在最前面) ===
fs = 200e6;
Start_loc_Base = 469200000;
Rcorr_Default = 0.35;

% --- 局部放大区域设置 (统一) ---
Zoom_Az_Lim = [165, 171];  % 方位角放大范围
Zoom_El_Lim = [40, 44];    % 仰角放大范围

% ==========================================
% 数据集 1: 传统方法 (Fixed Window / Loop)
% ==========================================
% 请确保这里的路径是你电脑上正确的文件路径
file_trad = '..\2024\2d\results\20240822165932_loop_result_yld_4.692e8_1.8e6_window_1024_256_ErrCalc_upsample_50.txt';
Th_Trad_Err = 0.7; % 传统方法的误差阈值 

% ==========================================
% 数据集 2: 本文方法 (Adaptive Window)
% ==========================================
% 请确保这里的路径是你电脑上正确的文件路径
file_prop = '..\2024\2d\results\result_yld_ADAPTIVE_upsample_469000000_472000000.txt';
Th_Prop_Err = 0.5; % 本文方法的误差阈值 

%% === 2. 数据读取与处理 ===
% --- 读取并筛选数据 1 (传统) ---
data_trad = process_data(file_trad, Start_loc_Base, fs, Th_Trad_Err, Rcorr_Default);
% --- 读取并筛选数据 2 (本文) ---
data_prop = process_data(file_prop, Start_loc_Base, fs, Th_Prop_Err, Rcorr_Default);

% --- 计算全局统一的时间轴范围 ---
min_t = min([min(data_trad.Time_us); min(data_prop.Time_us)]);
max_t = max([max(data_trad.Time_us); max(data_prop.Time_us)]);
Global_Time_Lim = [min_t, max_t];

% --- 计算全局统一的全景图范围 ---
min_az = min([min(data_trad.Azimuth); min(data_prop.Azimuth)]);
max_az = max([max(data_trad.Azimuth); max(data_prop.Azimuth)]);
min_el = min([min(data_trad.Elevation); min(data_prop.Elevation)]);
max_el = max([max(data_trad.Elevation); max(data_prop.Elevation)]);
Global_Overview_Az = [min_az-2, max_az+2];
Global_Overview_El = [min_el-2, max_el+2];

% =========================================================================
% === 2.1 统计误差计算与打印 (Statistics Calculation) ===
% =========================================================================
fprintf('\n================================================================================\n');
fprintf('                        统计分析结果 (Statistics Report)                        \n');
fprintf('================================================================================\n');

% 定义计算统计量的匿名函数
calc_stat_struct = @(x) struct('Mean', mean(x), 'Std', std(x), 'RMS', rms(x), 'Max', max(x));

% --- 计算传统方法统计量 ---
stat_trad_az = calc_stat_struct(data_trad.Err_Az);
stat_trad_el = calc_stat_struct(data_trad.Err_El);

% --- 计算本文方法统计量 ---
stat_prop_az = calc_stat_struct(data_prop.Err_Az);
stat_prop_el = calc_stat_struct(data_prop.Err_El);

% --- 打印方位角 (Azimuth) 对比 ---
fprintf('>> [Azimuth Error] (方位角误差统计)\n');
fprintf('%-20s | %-10s | %-10s | %-10s | %-10s\n', 'Method', 'Mean(°)', 'Std(°)', 'RMS(°)', 'Max(°)');
fprintf('---------------------|------------|------------|------------|------------\n');
fprintf('%-20s | %10.4f | %10.4f | %10.4f | %10.4f\n', 'Traditional', ...
    stat_trad_az.Mean, stat_trad_az.Std, stat_trad_az.RMS, stat_trad_az.Max);
fprintf('%-20s | %10.4f | %10.4f | %10.4f | %10.4f\n', 'Proposed', ...
    stat_prop_az.Mean, stat_prop_az.Std, stat_prop_az.RMS, stat_prop_az.Max);
% 计算提升比例 (以 RMS 为准)
if stat_trad_az.RMS > 0
    imp_az = (stat_trad_az.RMS - stat_prop_az.RMS) / stat_trad_az.RMS * 100;
    fprintf('** Azimuth RMS Improvement: %.2f%% **\n\n', imp_az);
else
    fprintf('\n');
end

% --- 打印仰角 (Elevation) 对比 ---
fprintf('>> [Elevation Error] (仰角误差统计)\n');
fprintf('%-20s | %-10s | %-10s | %-10s | %-10s\n', 'Method', 'Mean(°)', 'Std(°)', 'RMS(°)', 'Max(°)');
fprintf('---------------------|------------|------------|------------|------------\n');
fprintf('%-20s | %10.4f | %10.4f | %10.4f | %10.4f\n', 'Traditional', ...
    stat_trad_el.Mean, stat_trad_el.Std, stat_trad_el.RMS, stat_trad_el.Max);
fprintf('%-20s | %10.4f | %10.4f | %10.4f | %10.4f\n', 'Proposed', ...
    stat_prop_el.Mean, stat_prop_el.Std, stat_prop_el.RMS, stat_prop_el.Max);
% 计算提升比例
if stat_trad_el.RMS > 0
    imp_el = (stat_trad_el.RMS - stat_prop_el.RMS) / stat_trad_el.RMS * 100;
    fprintf('** Elevation RMS Improvement: %.2f%% **\n', imp_el);
else
    fprintf('\n');
end

% --- 打印有效点数对比 ---
fprintf('\n>> [Data Count] (有效点数)\n');
fprintf('Traditional: %d points\n', height(data_trad));
fprintf('Proposed:    %d points\n', height(data_prop));
fprintf('================================================================================\n\n');

%% === 3. 绘图 (4行 x 2列) ===
% 设置画布大小
fig_width = 20;   
fig_height = 26;  
f = figure('Units', 'centimeters', 'Position', [5, 2, fig_width, fig_height], 'Color', 'w');
% 字体设置
font_name = 'Microsoft YaHei'; font_size = 10; label_size = 11;
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
col_titles = {'传统方法', '本文方法'};

% --- 初始化字母计数器 ---
plot_count = 1;

% =========================================================================
% 第 1 行：二维全景图 (Overview) -> (a), (b)
% =========================================================================
ax1_1 = nexttile;
plot_overview(ax1_1, data_trad, Zoom_Az_Lim, Zoom_El_Lim, Global_Overview_Az, Global_Overview_El);
title(col_titles{1}, 'FontName', font_name, 'FontSize', label_size+2, 'FontWeight', 'bold');
apply_jgr_style(ax1_1, font_name, font_size, sprintf('(%s)', char(96 + plot_count))); plot_count = plot_count + 1;
ylabel('仰角 (°)', 'FontName', font_name, 'FontSize', label_size);
xlabel('方位角 (°)', 'FontName', font_name, 'FontSize', label_size);

ax1_2 = nexttile;
plot_overview(ax1_2, data_prop, Zoom_Az_Lim, Zoom_El_Lim, Global_Overview_Az, Global_Overview_El);
title(col_titles{2}, 'FontName', font_name, 'FontSize', label_size+2, 'FontWeight', 'bold');
apply_jgr_style(ax1_2, font_name, font_size, sprintf('(%s)', char(96 + plot_count))); plot_count = plot_count + 1;
xlabel('方位角 (°)', 'FontName', font_name, 'FontSize', label_size);

cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = '时间 (\mus)';
cb.Label.FontName = font_name; cb.Label.FontSize = label_size; colormap('jet');

% =========================================================================
% 第 2 行：局部放大图 (Zoomed View) -> (c), (d)
% =========================================================================
ax2_1 = nexttile;
plot_zoom(ax2_1, data_trad, Zoom_Az_Lim, Zoom_El_Lim);
apply_jgr_style(ax2_1, font_name, font_size, sprintf('(%s)', char(96 + plot_count))); plot_count = plot_count + 1;
ylabel('仰角 (°)', 'FontName', font_name, 'FontSize', label_size);
xlabel('方位角 (°)', 'FontName', font_name, 'FontSize', label_size);

ax2_2 = nexttile;
plot_zoom(ax2_2, data_prop, Zoom_Az_Lim, Zoom_El_Lim);
apply_jgr_style(ax2_2, font_name, font_size, sprintf('(%s)', char(96 + plot_count))); plot_count = plot_count + 1;
xlabel('方位角 (°)', 'FontName', font_name, 'FontSize', label_size);

% =========================================================================
% 第 3 行：方位角误差 (Azimuth Error) -> (e), (f)
% =========================================================================
ax3_1 = nexttile;
plot_error_time(ax3_1, data_trad, 'Err_Az', Global_Time_Lim, Th_Trad_Err);
apply_jgr_style(ax3_1, font_name, font_size, sprintf('(%s)', char(96 + plot_count))); plot_count = plot_count + 1;
ylabel('方位角误差 (°)', 'FontName', font_name, 'FontSize', label_size);
xlabel('时间 (\mus)', 'FontName', font_name, 'FontSize', label_size);

ax3_2 = nexttile;
plot_error_time(ax3_2, data_prop, 'Err_Az', Global_Time_Lim, Th_Trad_Err);
apply_jgr_style(ax3_2, font_name, font_size, sprintf('(%s)', char(96 + plot_count))); plot_count = plot_count + 1;
xlabel('时间 (\mus)', 'FontName', font_name, 'FontSize', label_size);

% =========================================================================
% 第 4 行：仰角误差 (Elevation Error) -> (g), (h)
% =========================================================================
ax4_1 = nexttile;
plot_error_time(ax4_1, data_trad, 'Err_El', Global_Time_Lim, Th_Trad_Err);
apply_jgr_style(ax4_1, font_name, font_size, sprintf('(%s)', char(96 + plot_count))); plot_count = plot_count + 1;
ylabel('仰角误差 (°)', 'FontName', font_name, 'FontSize', label_size);
xlabel('时间 (\mus)', 'FontName', font_name, 'FontSize', label_size);

ax4_2 = nexttile;
plot_error_time(ax4_2, data_prop, 'Err_El', Global_Time_Lim, Th_Trad_Err);
apply_jgr_style(ax4_2, font_name, font_size, sprintf('(%s)', char(96 + plot_count)));
xlabel('时间 (\mus)', 'FontName', font_name, 'FontSize', label_size);

linkaxes([ax3_1, ax3_2, ax4_1, ax4_2], 'x');

%% === 辅助函数定义 ===

% 1. 数据处理函数
function data = process_data(fname, base_loc, fs, err_th, rcorr_def)
    if ~isfile(fname), error(['文件不存在: ' fname]); end
    opts = detectImportOptions(fname); opts.VariableNamingRule = 'preserve';
    T = readtable(fname, opts);
    
    % 基础地理筛选
    valid_geo = T.Start_loc > base_loc & T.Start_loc < (base_loc + 1800000) & ...
                T.Elevation < 85 & T.Azimuth < 250 & abs(T.t123) < 1.0;
            
    % 坏点区域剔除
    bad_region = false(height(T), 1);
    if ismember('Azimuth', T.Properties.VariableNames)
        bad_region = (T.Azimuth > 160) & (T.Azimuth < 250) & (T.Elevation > 0) & (T.Elevation < 35);
    end
    
    % Rcorr 筛选
    mask_rcorr = false(height(T), 1);
    if ismember('Win_Len', T.Properties.VariableNames)
        idx_512 = (T.Win_Len == 512);   mask_rcorr(idx_512) = T.Rcorr(idx_512) > 0.1;
        idx_1024 = (T.Win_Len == 1024); mask_rcorr(idx_1024) = T.Rcorr(idx_1024) > 0.1; 
        idx_2048 = (T.Win_Len == 2048); mask_rcorr(idx_2048) = T.Rcorr(idx_2048) > 0.1;
        idx_4096 = (T.Win_Len == 4096); mask_rcorr(idx_4096) = T.Rcorr(idx_4096) > 0.1;
        other = ~ismember(T.Win_Len, [512, 1024, 2048, 4096]);
        mask_rcorr(other) = T.Rcorr(other) > rcorr_def;
    else
        mask_rcorr = T.Rcorr > rcorr_def;
    end
    
    % 误差筛选
    mask_error = (T.Err_Az < err_th) & (T.Err_El < err_th);
    
    final_idx = valid_geo & (~bad_region) & mask_rcorr & mask_error;
    data = T(final_idx, :);
    data.Time_us = (data.Start_loc - base_loc) / fs * 1e6;
end

% 2. 绘图函数 - Overview
function plot_overview(ax, data, zoom_az, zoom_el, limit_az, limit_el)
    axes(ax); hold on; axis equal; grid on;
    errorbar(data.Azimuth, data.Elevation, data.Err_El, data.Err_El, data.Err_Az, data.Err_Az, ...
        '.', 'Color', [0.8 0.8 0.8], 'CapSize', 0, 'LineWidth', 0.1); 
    scatter(data.Azimuth, data.Elevation, 4, data.Time_us, 'filled');
    rectangle('Position', [zoom_az(1), zoom_el(1), diff(zoom_az), diff(zoom_el)], ...
              'EdgeColor', 'r', 'LineWidth', 1.2, 'LineStyle', '-');
    xlim(limit_az); ylim(limit_el);
    caxis([min(data.Time_us), max(data.Time_us)]);
end

% 3. 绘图函数 - Zoom
function plot_zoom(ax, data, zoom_az, zoom_el)
    axes(ax); hold on; axis equal; grid on;
    scatter(data.Azimuth, data.Elevation, 15, data.Time_us, 'filled');
    eb = errorbar(data.Azimuth, data.Elevation, data.Err_El, data.Err_El, data.Err_Az, data.Err_Az, ...
        '.', 'Color', [0.4 0.4 0.4], 'CapSize', 0, 'LineWidth', 1.0);
    eb.Annotation.LegendInformation.IconDisplayStyle = 'off';
    xlim(zoom_az); ylim(zoom_el);
    set(ax, 'XColor', 'r', 'YColor', 'r', 'LineWidth', 1.2);
end

% 4. 修改后的误差随时间变化绘图函数 (Mean和RMS移至右侧)
function plot_error_time(ax, data, field, t_lim, y_lim_val)
    axes(ax); hold on; grid on;
    scatter(data.Time_us, data.(field), 8, data.Time_us, 'filled', 'MarkerFaceAlpha', 0.6);
    xlim(t_lim); 
    ylim([0, y_lim_val]);
    
end

% 5. 样式应用 (保持不变，但确保调用时传入的是 (a), (b) 等)
function apply_jgr_style(ax, fname, fsize, label_text)
    set(ax, 'FontName', fname, 'FontSize', fsize, 'LineWidth', 1.0, ...
        'Box', 'on', 'TickDir', 'in', 'GridAlpha', 0.3, 'GridLineStyle', ':');
    % 将子图标签放置在左上角
    text(ax, 0.03, 0.93, label_text, 'Units', 'normalized', ...
        'FontName', fname, 'FontSize', fsize+1, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
end