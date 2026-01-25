clear; clc; close all;

%% ================== 1. 用户参数设置 ==================
filename = 'results\20230618125747.5480_400000000_99999999_1024_256_8_gage-20230306.txt';
r_loction_yld = 400955647;       
signal_length = 3.7e5;         
sampling_interval_ns = 5;      
ns_to_us = 1e-3; 

% --- 局部放大区域控制 ---
Zoom_Start_us = 990;    
Zoom_End_us   = 1190;   

% --- 分支分类阈值 ---
Azimuth_Split_Threshold = 300; 

% --- 【关键】背景色块的时间分辨率 ---
% 越小越精细，越大越平滑。建议设为 1~5 us，取决于先导步进的频率
Time_Bin_Size = 1.3; 

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
vhf_time = (filteredTable.Start_loc - r_loction_yld) * time_conversion_factor;
vhf_el   = filteredTable.Elevation;
vhf_az   = filteredTable.Azimuth; 

%% ================== 3. 分支分离 ==================
zoom_idx = (vhf_time >= Zoom_Start_us) & (vhf_time <= Zoom_End_us);
z_time = vhf_time(zoom_idx);
z_el   = vhf_el(zoom_idx);
z_az   = vhf_az(zoom_idx);

mask_b1 = z_az < Azimuth_Split_Threshold; % 分支 A (红色)
b1_time = z_time(mask_b1);
b1_el   = z_el(mask_b1);
b1_az   = z_az(mask_b1);

mask_b2 = z_az >= Azimuth_Split_Threshold; % 分支 B (蓝色)
b2_time = z_time(mask_b2);
b2_el   = z_el(mask_b2);
b2_az   = z_az(mask_b2);

%% ================== 4. 绘图 ==================
figure('Units', 'centimeters', 'Position', [5, 5, 24, 10], 'Color', 'w');
t = tiledlayout(1, 3, 'TileSpacing', 'tight', 'Padding', 'compact');

font_name = 'Arial';
font_size = 10;
sz_fg = 25;      

% === 子图 1: 全局空间视角 (无分割线) ===
ax1 = nexttile; 
hold(ax1, 'on');

% 1.1 绘制全背景 (浅灰)
scatter(ax1, vhf_az, vhf_el, 10, [0.85 0.85 0.85], 'filled');

% 1.2 绘制两个分支 (红/蓝)
scatter(ax1, b1_az, b1_el, sz_fg, 'r', 'filled', 'MarkerFaceAlpha', 0.8);
scatter(ax1, b2_az, b2_el, sz_fg, 'b', 'filled', 'MarkerFaceAlpha', 0.8);

% 1.3 红色包围框
margin = 0.5; 
if ~isempty(z_az)
    rect_x = min(z_az) - margin; 
    rect_y = min(z_el) - margin;
    rect_w = (max(z_az) + margin) - rect_x;
    rect_h = (max(z_el) + margin) - rect_y;
    rectangle(ax1, 'Position', [rect_x, rect_y, rect_w, rect_h], ...
        'EdgeColor', 'k', 'LineWidth', 1.0, 'LineStyle', '--');
end

ylabel(ax1, 'Elevation (°)', 'FontName', font_name, 'FontSize', font_size);
xlabel(ax1, 'Azimuth (°)', 'FontName', font_name, 'FontSize', font_size);
title(ax1, '(a) 2D Spatial Overview', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
apply_jgr_style(ax1, font_name, font_size);
axis(ax1, 'equal'); axis(ax1, 'tight');
ylim(ax1, [min(z_el)-5, max(z_el)+5]); 

% === 子图 2: 时空交替演化 (带背景色块) ===
ax2 = nexttile([1, 2]); 
hold(ax2, 'on');

% --- Step A: 绘制背景互斥色块 (核心步骤) ---
% 定义Y轴覆盖范围 (铺满整个Y轴)
if ~isempty(z_el)
    y_min_fill = min(z_el) - 2;
    y_max_fill = max(z_el) + 2;
else
    y_min_fill = 0; y_max_fill = 90;
end

% 定义时间网格
time_grid = Zoom_Start_us : Time_Bin_Size : Zoom_End_us;

% 循环检测每个时间窗口
for i = 1:length(time_grid)-1
    t_s = time_grid(i);
    t_e = time_grid(i+1);
    
    % 检测该时间段内，两个分支是否有辐射源点
    has_b1 = any(b1_time >= t_s & b1_time < t_e);
    has_b2 = any(b2_time >= t_s & b2_time < t_e);
    
    fill_color = [];
    
    if has_b1 && ~has_b2
        % 只有 A 活动 -> 淡红色背景
        fill_color = [1, 0.92, 0.92]; 
    elseif has_b2 && ~has_b1
        % 只有 B 活动 -> 淡蓝色背景
        fill_color = [0.92, 0.92, 1]; 
    elseif has_b1 && has_b2
        % 同时活动 (非常少见，如果有可以用紫色或灰色)
        % fill_color = [0.95, 0.9, 0.95]; 
    end
    
    if ~isempty(fill_color)
        % 绘制色块 (不画边框 EdgeColor='none')
        fill(ax2, [t_s, t_e, t_e, t_s], [y_min_fill, y_min_fill, y_max_fill, y_max_fill], ...
            fill_color, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

% --- Step B: 绘制连线 (保留线条以显示阶梯感) ---
% 使用 stairs (阶梯图) 或 plot (折线图)。这里用 plot 配合点更自然
if ~isempty(b1_time)
    plot(ax2, b1_time, b1_el, '-', 'Color', [0.8 0 0, 0.3], 'LineWidth', 1.0); 
end
if ~isempty(b2_time)
    plot(ax2, b2_time, b2_el, '-', 'Color', [0 0 0.8, 0.3], 'LineWidth', 1.0); 
end

% --- Step C: 绘制散点 (在最上层) ---
s1 = scatter(ax2, b1_time, b1_el, sz_fg, 'r', 'filled', 'MarkerFaceAlpha', 1);
s2 = scatter(ax2, b2_time, b2_el, sz_fg, 'b', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 1); % 蓝色分支加黑边区分形状

% 样式调整
xlim(ax2, [Zoom_Start_us, Zoom_End_us]);
ylim(ax2, [y_min_fill, y_max_fill]);

ylabel(ax2, 'Elevation (°)', 'FontName', font_name, 'FontSize', font_size);
xlabel(ax2, 'Time (\mus)', 'FontName', font_name, 'FontSize', font_size);
title(ax2, '(b) Stepwise & Alternating Propagation', 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');

% 图例
legend([s1, s2], {'Branch 1 ', 'Branch 2 '}, ...
    'Location', 'northwest', 'Box', 'on', 'FontSize', 9);

apply_jgr_style(ax2, font_name, font_size);
set(ax2, 'Layer', 'top'); % 确保坐标轴刻度在色块之上
set(ax2, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1.2); 

fprintf('绘图完成。\n');

%% ================== 辅助函数 ==================
function apply_jgr_style(ax, fname, fsize)
    set(ax, 'FontName', fname, 'FontSize', fsize);
    set(ax, 'LineWidth', 1.0);      
    set(ax, 'TickDir', 'in'); % 刻度朝外，避免被色块遮挡      
    set(ax, 'Box', 'on');           
    set(ax, 'GridAlpha', 0.3);      
    set(ax, 'GridLineStyle', ':'); 
    grid(ax, 'on');
end