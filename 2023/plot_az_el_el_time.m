clear; clc; close all;

%% ================== 1. 用户参数设置 ==================
filename = 'results\20230618125747.5480_400000000_99999999_1024_256_8_gage-20230306.txt';
r_loction_yld = 400955647;       
signal_length = 3.7e5;         
sampling_interval_ns = 5;      
ns_to_us = 1e-3; 

% --- 【关键】局部放大区域控制 (单位: us) ---
Zoom_Start_us = 990;   % 放大起始时间
Zoom_End_us   = 1190;   % 放大结束时间

%% ================== 2. 数据加载与预处理 ==================
if ~isfile(filename), error('文件不存在: %s', filename); end
fprintf('正在加载数据...\n');
result1 = readtable(filename);

time_conversion_factor = sampling_interval_ns * ns_to_us; 
filter_loc_min = r_loction_yld;
filter_loc_max = r_loction_yld + signal_length;

logicalIndex =  abs(result1.t123) < 1  & ...
                abs(result1.Rcorrn) > 0.3 & ...
                result1.Start_loc < filter_loc_max & ...
                result1.Start_loc > filter_loc_min & ...
                result1.Elevation < 80 & ...
                result1.Elevation > 0 & ...
                result1.Azimuth > 255 & ...
                result1.Azimuth < 330;

filteredTable = result1(logicalIndex, :);

% 计算物理量
vhf_time = (filteredTable.Start_loc - r_loction_yld) * time_conversion_factor;
vhf_el   = filteredTable.Elevation;
vhf_az   = filteredTable.Azimuth; 
vhf_color = vhf_time; 

%% ================== 3. JGR 风格绘图 ==================
% 设置画布 (宽一些，比如 20cm，以容纳左右布局)
figure('Units', 'centimeters', 'Position', [5, 5, 20, 9], 'Color', 'w');

% 【关键布局修改】：1行3列
t = tiledlayout(1, 3, 'TileSpacing', 'tight', 'Padding', 'compact');

font_name = 'Arial';
font_size = 10;
point_size = 5;

% === 子图 1: 空间形态 (Spatial View) - 占 1/3 宽度 ===
ax1 = nexttile; 
scatter(ax1, vhf_az, vhf_el, point_size, vhf_color, 'filled', 'MarkerFaceAlpha', 0.8);
hold on;

% --- 绘制红色空间包围框 ---
in_zoom_mask = (vhf_time >= Zoom_Start_us) & (vhf_time <= Zoom_End_us);

if any(in_zoom_mask)
    % 获取放大区域的空间边界
    az_zoom = vhf_az(in_zoom_mask);
    el_zoom = vhf_el(in_zoom_mask);
    
    margin = 0.5; 
    min_az = min(az_zoom) - margin; rect_w = (max(az_zoom) + margin) - min_az;
    min_el = min(el_zoom) - margin; rect_h = (max(el_zoom) + margin) - min_el;
    
    rectangle(ax1, 'Position', [min_az, min_el, rect_w, rect_h], ...
        'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '-');
end

% 样式
ylabel(ax1, 'Elevation (°)', 'FontName', font_name, 'FontSize', font_size);
xlabel(ax1, 'Azimuth (°)', 'FontName', font_name, 'FontSize', font_size);
title(ax1, '(a) 2D Spatial View', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
apply_jgr_style(ax1, font_name, font_size);
axis(ax1, 'equal'); % 保持空间比例一致 (可选，如果图形太扁可注释掉)
axis(ax1, 'tight');

% === 子图 2: 时间演化 (Time View) - 占 2/3 宽度 ===
% 【关键】：跨越 2 列 (1行, 2列)
ax2 = nexttile([1, 2]); 
scatter(ax2, vhf_time, vhf_el, point_size+10, vhf_color, 'filled', 'MarkerFaceAlpha', 0.8);

% 限制范围
xlim(ax2, [Zoom_Start_us, Zoom_End_us]);
xticks(ax2, Zoom_Start_us:5:Zoom_End_us);
if any(in_zoom_mask)
    y_data_view = vhf_el(in_zoom_mask);
    margin_y = 1; 
    ylim(ax2, [min(y_data_view)-margin_y, max(y_data_view)+margin_y]);
end

% 样式
% 右图的 Y 轴也是 Elevation，为了节省空间，可以选择保留或隐藏 Y 轴标签
ylabel(ax2, 'Elevation (°)', 'FontName', font_name, 'FontSize', font_size);
xlabel(ax2, 'Time (\mus)', 'FontName', font_name, 'FontSize', font_size);
title(ax2, '(b) Zoomed Temporal Evolution', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
apply_jgr_style(ax2, font_name, font_size);

% 给右图加框
set(ax2, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1.2); 

% === 公共 Colorbar ===
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Time (\mus)';
cb.Label.FontName = font_name;
cb.Label.FontSize = font_size;
colormap('jet');

fprintf('绘图完成。\n');

%% ================== 辅助函数 ==================
function apply_jgr_style(ax, fname, fsize)
    set(ax, 'FontName', fname, 'FontSize', fsize);
    set(ax, 'LineWidth', 1.0);      
    set(ax, 'TickDir', 'in');       
    set(ax, 'Box', 'on');           
    set(ax, 'GridAlpha', 0.2);      
    set(ax, 'GridLineStyle', '--'); 
    grid(ax, 'on');
end