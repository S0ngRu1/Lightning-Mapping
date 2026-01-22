clear; clc; close all;

%% === 1. 全局参数与文件设置 ===
fs = 200e6;
Start_loc_Base = 469200000;
Rcorr_Default = 0.35;

% --- 局部放大区域设置 (统一) ---
Zoom_Az_Lim = [167, 177];  % 方位角放大范围
Zoom_El_Lim = [39.5, 45.5];    % 仰角放大范围

% ==========================================
% 数据集 1: 传统方法 (Fixed Window / Loop)
% ==========================================
file_trad = '..\2024\2d\results\20240822165932_loop_result_yld_4.692e8_1.8e6_window_1024_256_ErrCalc_upsample_50.txt';
Th_Trad_Err = 0.7; % 传统方法的误差阈值 

% ==========================================
% 数据集 2: 本文方法 (Adaptive Window)
% ==========================================
file_prop = '..\2024\2d\results\result_yld_3.65e8_5.6e8_window_ADAPTIVE_1e4_factor3_with_error.txt';
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

%% === 3. 绘图 (4行 x 2列) ===
% 设置画布大小 (A4 比例拉长)
fig_width = 20;   % 宽 (cm)
fig_height = 24;  % 高 (cm)
f = figure('Units', 'centimeters', 'Position', [5, 2, fig_width, fig_height], 'Color', 'w');

% 字体设置
font_name = 'Arial';
font_size = 10;
label_size = 11;

% 布局
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 定义列标题
col_titles = {'Traditional Method', 'Proposed Method'};

% =========================================================================
% 第 1 行：二维全景图 (Overview)
% =========================================================================
% --- Col 1: Traditional ---
ax1_1 = nexttile;
plot_overview(ax1_1, data_trad, Zoom_Az_Lim, Zoom_El_Lim, Global_Overview_Az, Global_Overview_El);
title(col_titles{1}, 'FontName', font_name, 'FontSize', label_size+2, 'FontWeight', 'bold');
apply_jgr_style(ax1_1, font_name, font_size, '(a1)');
ylabel('Elevation (°)', 'FontName', font_name, 'FontSize', label_size);

% --- Col 2: Proposed ---
ax1_2 = nexttile;
plot_overview(ax1_2, data_prop, Zoom_Az_Lim, Zoom_El_Lim, Global_Overview_Az, Global_Overview_El);
title(col_titles{2}, 'FontName', font_name, 'FontSize', label_size+2, 'FontWeight', 'bold');
apply_jgr_style(ax1_2, font_name, font_size, '(a2)');
% 右侧添加 Colorbar
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Time (\mus)';
cb.Label.FontName = font_name; cb.Label.FontSize = label_size;
colormap('jet');

% =========================================================================
% 第 2 行：局部放大图 (Zoomed View)
% =========================================================================
% --- Col 1: Traditional ---
ax2_1 = nexttile;
plot_zoom(ax2_1, data_trad, Zoom_Az_Lim, Zoom_El_Lim);
apply_jgr_style(ax2_1, font_name, font_size, '(b1)');
ylabel('Elevation (°)', 'FontName', font_name, 'FontSize', label_size);
% 【修改】添加横坐标标题
xlabel('Azimuth (°)', 'FontName', font_name, 'FontSize', label_size);

% --- Col 2: Proposed ---
ax2_2 = nexttile;
plot_zoom(ax2_2, data_prop, Zoom_Az_Lim, Zoom_El_Lim);
apply_jgr_style(ax2_2, font_name, font_size, '(b2)');
% 【修改】添加横坐标标题
xlabel('Azimuth (°)', 'FontName', font_name, 'FontSize', label_size);

% =========================================================================
% 第 3 行：方位角误差 (Azimuth Error vs Time)
% =========================================================================
% --- Col 1: Traditional ---
ax3_1 = nexttile;
plot_error_time(ax3_1, data_trad, 'Err_Az', Global_Time_Lim, Th_Trad_Err);
apply_jgr_style(ax3_1, font_name, font_size, '(c1)');
ylabel('Azimuth Error (°)', 'FontName', font_name, 'FontSize', label_size);
xticklabels([]); % 中间行不显示X轴刻度

% --- Col 2: Proposed ---
ax3_2 = nexttile;
plot_error_time(ax3_2, data_prop, 'Err_Az', Global_Time_Lim, Th_Prop_Err); % 注意阈值不同
apply_jgr_style(ax3_2, font_name, font_size, '(c2)');
xticklabels([]);

% =========================================================================
% 第 4 行：仰角误差 (Elevation Error vs Time)
% =========================================================================
% --- Col 1: Traditional ---
ax4_1 = nexttile;
plot_error_time(ax4_1, data_trad, 'Err_El', Global_Time_Lim, Th_Trad_Err);
apply_jgr_style(ax4_1, font_name, font_size, '(d1)');
ylabel('Elevation Error (°)', 'FontName', font_name, 'FontSize', label_size);
xlabel('Time (\mus)', 'FontName', font_name, 'FontSize', label_size);

% --- Col 2: Proposed ---
ax4_2 = nexttile;
plot_error_time(ax4_2, data_prop, 'Err_El', Global_Time_Lim, Th_Prop_Err);
apply_jgr_style(ax4_2, font_name, font_size, '(d2)');
xlabel('Time (\mus)', 'FontName', font_name, 'FontSize', label_size);

% 链接时间轴，方便缩放
linkaxes([ax3_1, ax3_2, ax4_1, ax4_2], 'x');


%% === 辅助函数定义 ===

% 1. 数据处理函数
function data = process_data(fname, base_loc, fs, err_th, rcorr_def)
    if ~isfile(fname), error(['文件不存在: ' fname]); end
    opts = detectImportOptions(fname); opts.VariableNamingRule = 'preserve';
    T = readtable(fname, opts);
    
    % 基础地理筛选
    valid_geo = T.Start_loc > base_loc & T.Start_loc < (base_loc + 1800000) & ...
                T.Elevation < 85 & T.Azimuth < 250 & abs(T.t123) < 10.0;
            
    % 坏点区域剔除
    bad_region = false(height(T), 1);
    if ismember('Azimuth', T.Properties.VariableNames)
        bad_region = (T.Azimuth > 160) & (T.Azimuth < 250) & (T.Elevation > 0) & (T.Elevation < 35);
    end
    
    % Rcorr 筛选
    mask_rcorr = false(height(T), 1);
    if ismember('Win_len', T.Properties.VariableNames)
        idx_512 = (T.Win_len == 512);   mask_rcorr(idx_512) = T.Rcorr(idx_512) > 0.3;
        idx_1024 = (T.Win_len == 1024); mask_rcorr(idx_1024) = T.Rcorr(idx_1024) > 0.25; 
        idx_2048 = (T.Win_len == 2048); mask_rcorr(idx_2048) = T.Rcorr(idx_2048) > 0.1;
        idx_4096 = (T.Win_len == 4096); mask_rcorr(idx_4096) = T.Rcorr(idx_4096) > 0.1;
        other = ~ismember(T.Win_len, [512, 1024, 2048, 4096]);
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
    % 误差棒
    errorbar(data.Azimuth, data.Elevation, data.Err_El, data.Err_El, data.Err_Az, data.Err_Az, ...
        '.', 'Color', [0.8 0.8 0.8], 'CapSize', 0, 'LineWidth', 0.1); 
    % 散点
    scatter(data.Azimuth, data.Elevation, 4, data.Time_us, 'filled');
    % 红框
    rectangle('Position', [zoom_az(1), zoom_el(1), diff(zoom_az), diff(zoom_el)], ...
              'EdgeColor', 'r', 'LineWidth', 1.2, 'LineStyle', '-');
    xlim(limit_az); ylim(limit_el);
    caxis([min(data.Time_us), max(data.Time_us)]); % 确保颜色一致
end

% 3. 绘图函数 - Zoom
function plot_zoom(ax, data, zoom_az, zoom_el)
    axes(ax); hold on; axis equal; grid on;
    % 误差棒
    eb = errorbar(data.Azimuth, data.Elevation, data.Err_El, data.Err_El, data.Err_Az, data.Err_Az, ...
        '.', 'Color', [0.6 0.6 0.6], 'CapSize', 0, 'LineWidth', 0.8);
    eb.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % 散点
    scatter(data.Azimuth, data.Elevation, 15, data.Time_us, 'filled');
    xlim(zoom_az); ylim(zoom_el);
    set(ax, 'XColor', 'r', 'YColor', 'r', 'LineWidth', 1.2); % 红色边框
end

% 4. 绘图函数 - Error vs Time
function plot_error_time(ax, data, field, t_lim, y_lim_val)
    axes(ax); hold on; grid on;
    scatter(data.Time_us, data.(field), 8, data.Time_us, 'filled', 'MarkerFaceAlpha', 0.6);
    xlim(t_lim); 
    ylim([0, y_lim_val]);
end

% 5. 样式应用
function apply_jgr_style(ax, fname, fsize, label_text)
    set(ax, 'FontName', fname, 'FontSize', fsize, 'LineWidth', 1.0, ...
        'Box', 'on', 'TickDir', 'in', 'GridAlpha', 0.3, 'GridLineStyle', ':');
    text(ax, 0.03, 0.93, label_text, 'Units', 'normalized', ...
        'FontName', fname, 'FontSize', fsize+1, 'FontWeight', 'bold');
end