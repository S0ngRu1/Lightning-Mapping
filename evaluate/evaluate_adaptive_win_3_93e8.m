clear; clc; close all;

%% === 1. 参数设置 ===
% --- 文件路径 ---
file_col1 = '..\2024\2d\results\20240822165932_result_yld_3.65e8_5e8_window_512_128_阈值4倍标准差_去零飘_30_80_hann_with_error.txt';
file_col2 = '..\2024\2d\results\20240822165932_result_yld_3.65e8_5e8_window_4096_1024_阈值4倍标准差_去零飘_30_80_hann.txt';
file_col3 = '..\2024\2d\results\result_yld_3.65e8_5.6e8_window_ADAPTIVE_1e4_factor3.txt';

% --- 通用参数 ---
fs = 200e6;              % 采样率
Start_loc_Base = 3.93e8;  % 基准采样点位置
T_us_max = 35000;        % 最大显示时间 (us)

% --- 【关键】局部放大区域设置 ---
Zoom_Az_Lim = [140, 148];  % 方位角放大范围
Zoom_El_Lim = [30, 42];    % 仰角放大范围

% --- 全局视图范围 ---
Global_Az_Lim = [125, 170];
Global_El_Lim = [5, 45];

%% === 2. 数据读取与预处理 ===
% Col 1: 512窗口
data_1 = read_and_filter(file_col1, Start_loc_Base, 1, 0.6, T_us_max);
% Col 2: 4096窗口
data_2 = read_and_filter(file_col2, Start_loc_Base, 1, 0.1, T_us_max);
% Col 3: Adaptive
data_3 = read_and_filter(file_col3, Start_loc_Base, 1.0, 0.1, T_us_max);

%% === 3. 绘图 (5行 x 3列 虚拟网格) ===
% 设置画布大小
f = figure('Units', 'pixels', 'Color', 'w', 'Position', [50, 50, 1200, 700]);

% --- 【兼容性修改】使用 5x3 网格来实现 3:2 的高度比例 ---
t = tiledlayout(5, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% 注意：删除了 t.RowHeight 属性设置，以兼容 R2021a 以下版本

% 绘图参数
pt_size_overview = 2;   % 全景图点大小     
pt_size_zoom = 10;      % 放大图点大小
alpha_val = 0.8;    

% 定义列标题
col_titles = {'(win:512)', '(win:4096)', '(Adaptive win)'};

% =========================================================================
% 第 1 行：整体对比 (Overview) [占用前 3 行网格]
% =========================================================================

% --- 1.1 Col 1 Overview ---
% nexttile(索引, [跨行数, 跨列数])
ax1 = nexttile(1, [3, 1]); 
plot_overview(ax1, data_1, Zoom_Az_Lim, Zoom_El_Lim, Global_Az_Lim, Global_El_Lim, pt_size_overview, alpha_val);
title(['Overview ' col_titles{1}], 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Elevation (°)'); 
text(ax1, 0.02, 0.95, '(a1)', 'Units', 'normalized', 'FontWeight', 'bold');

% --- 1.2 Col 2 Overview ---
ax2 = nexttile(2, [3, 1]);
plot_overview(ax2, data_2, Zoom_Az_Lim, Zoom_El_Lim, Global_Az_Lim, Global_El_Lim, pt_size_overview, alpha_val);
title(['Overview ' col_titles{2}], 'FontSize', 11, 'FontWeight', 'bold');
text(ax2, 0.02, 0.95, '(a2)', 'Units', 'normalized', 'FontWeight', 'bold');

% --- 1.3 Col 3 Overview ---
ax3 = nexttile(3, [3, 1]);
plot_overview(ax3, data_3, Zoom_Az_Lim, Zoom_El_Lim, Global_Az_Lim, Global_El_Lim, pt_size_overview, alpha_val);
title(['Overview ' col_titles{3}], 'FontSize', 11, 'FontWeight', 'bold');
text(ax3, 0.02, 0.95, '(a3)', 'Units', 'normalized', 'FontWeight', 'bold');

% =========================================================================
% 第 2 行：局部放大对比 (Zoomed) [占用后 2 行网格]
% =========================================================================
% 计算起始索引：前3行共 3*3=9 个格子，所以第4行从索引 10 开始

% --- 2.1 Col 1 Zoom ---
ax4 = nexttile(10, [2, 1]);
plot_zoom(ax4, data_1, Zoom_Az_Lim, Zoom_El_Lim, pt_size_zoom, alpha_val);
ylabel('Elevation (°)'); xlabel('Azimuth (°)');
text(ax4, 0.02, 0.95, '(b1) Zoomed', 'Units', 'normalized', 'FontWeight', 'bold', 'Color', 'r');

% --- 2.2 Col 2 Zoom ---
ax5 = nexttile(11, [2, 1]);
plot_zoom(ax5, data_2, Zoom_Az_Lim, Zoom_El_Lim, pt_size_zoom, alpha_val);
xlabel('Azimuth (°)');
text(ax5, 0.02, 0.95, '(b2) Zoomed', 'Units', 'normalized', 'FontWeight', 'bold', 'Color', 'r');

% --- 2.3 Col 3 Zoom ---
ax6 = nexttile(12, [2, 1]);
plot_zoom(ax6, data_3, Zoom_Az_Lim, Zoom_El_Lim, pt_size_zoom, alpha_val);
xlabel('Azimuth (°)');
text(ax6, 0.02, 0.95, '(b3) Zoomed', 'Units', 'normalized', 'FontWeight', 'bold', 'Color', 'r');

% =========================================================================
% 全局设置
% =========================================================================
% 统一 Colorbar
cb = colorbar;
cb.Layout.Tile = 'east';
ylabel(cb, 'Time (\mus)', 'FontSize', 11);
colormap(jet);
caxis([0, T_us_max]); 

% 统一坐标轴风格
all_axes = findall(f, 'type', 'axes');
set(all_axes, 'FontSize', 10, 'LineWidth', 1.0, 'Box', 'on', ...
    'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15], 'TickDir', 'in');

fprintf('绘图完成。\n点数统计:\n Win512: %d\n Win4096: %d\n Adaptive: %d\n', ...
    height(data_1), height(data_2), height(data_3));

%% === 4. 辅助绘图函数 ===

% --- 绘制全景图 (带红框) ---
function plot_overview(ax, data, zoom_az, zoom_el, limit_az, limit_el, sz, alp)
    axes(ax); hold on; grid on;
    if ~isempty(data)
        scatter(data.Azimuth, data.Elevation, sz, data.Time_us, 'filled', 'MarkerFaceAlpha', alp);
    end
    % 画红框
    rectangle('Position', [zoom_az(1), zoom_el(1), diff(zoom_az), diff(zoom_el)], ...
              'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '-');
    xlim(limit_az); 
    ylim(limit_el);
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);
end

% --- 绘制局部放大图 (红色边框) ---
function plot_zoom(ax, data, zoom_az, zoom_el, sz, alp)
    axes(ax); hold on; grid on;
    if ~isempty(data)
        scatter(data.Azimuth, data.Elevation, sz, 'r', 'filled', 'MarkerFaceAlpha', alp);
    end
    xlim(zoom_az); 
    ylim(zoom_el);
    
    % 绘制红色矩形框
    rectangle('Position', [zoom_az(1), zoom_el(1), diff(zoom_az), diff(zoom_el)], ...
              'EdgeColor', 'r', 'LineWidth', 1, 'LineStyle', '-');
          
    % 设置网格线样式
    set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.5);
end

% --- 数据读取与筛选函数 (保持不变) ---
function filteredT = read_and_filter(fname, base_loc, t123, rcorr_default, T_us_max)
    if ~isfile(fname)
        warning(['文件不存在: ' fname]);
        filteredT = table();
        return;
    end
    
    opts = detectImportOptions(fname);
    opts.VariableNamingRule = 'preserve';
    T = readtable(fname, opts);
    
    % 1. 基本范围筛选
    max_samples = T_us_max * 200; 
    idx = T.Start_loc > base_loc & T.Start_loc < (base_loc + max_samples) ;
    if ismember('Azimuth', T.Properties.VariableNames) && ismember('Elevation', T.Properties.VariableNames)
        bad_region = (T.Azimuth > 160) & (T.Azimuth < 250) & ...
                     (T.Elevation > 0) & (T.Elevation < 50);
        idx = idx & (~bad_region);
    end
    
    % 2. 动态 Rcorr 筛选
    if ismember('Rcorr', T.Properties.VariableNames)
        if ismember('Win_Len', T.Properties.VariableNames)
            mask_rcorr = false(height(T), 1);
            
            idx_512 = (T.Win_Len == 512);
            mask_rcorr(idx_512) = T.Rcorr(idx_512) > 0.6;
            
            idx_1024 = (T.Win_Len == 1024);
            mask_rcorr(idx_1024) = T.Rcorr(idx_1024) > 0.4;
            
            idx_2048 = (T.Win_Len == 2048);
            mask_rcorr(idx_2048) = T.Rcorr(idx_2048) > 0.2;
            
            idx_4096 = (T.Win_Len == 4096);
            mask_rcorr(idx_4096) = T.Rcorr(idx_4096) > 0.1;
            
            other_lens = ~ismember(T.Win_Len, [512, 1024, 2048, 4096]);
            if any(other_lens)
                mask_rcorr(other_lens) = T.Rcorr(other_lens) > rcorr_default;
            end
            
            idx = idx & mask_rcorr;
        else
            idx = idx & (T.Rcorr > rcorr_default);
        end
    end
    
    % 3. 其他质量控制
    if ismember('t123', T.Properties.VariableNames)
        idx = idx & (abs(T.t123) < t123); 
    end
    if ismember('Elevation', T.Properties.VariableNames)
        idx = idx & (T.Elevation < 50);
    end
    if ismember('Azimuth', T.Properties.VariableNames)
        idx = idx & (T.Azimuth < 250);
    end
    
    filteredT = T(idx, :);
    filteredT.Time_us = (filteredT.Start_loc - base_loc) / 200e6 * 1e6;
end