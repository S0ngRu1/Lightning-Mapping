clear; clc; close all;

%% ==================== 1. 参数设置 ====================
filename = 'results\result_yld_3.65e8_5.6e8_window_ADAPTIVE_1e4_factor3_with_error.txt';
fs = 200e6;              
Error_Threshold = 1;   
Rcorr_Default = 0.3;     

% --- 时间/位置区间 ---
Event1_Range = [3.816e8, 3.88e8]; % 正先导 (Post-RS)
Event2_Range = [3.66e8, 3.72e8];  % 负先导 (Pre-RS)

%% ==================== 2. 数据读取与筛选 ====================
if ~isfile(filename), error('文件不存在: %s', filename); end
opts = detectImportOptions(filename); opts.VariableNamingRule = 'preserve';
T = readtable(filename, opts);

valid_geo = T.Start_loc > 3.6e8 & T.Start_loc < 4.0e8 & ...
            T.Elevation < 85 & T.Azimuth < 250 & abs(T.t123) < 1;

bad_region = false(height(T), 1);
if ismember('Azimuth', T.Properties.VariableNames)
    bad_region = (T.Azimuth > 160) & (T.Azimuth < 250) & (T.Elevation > 0) & (T.Elevation < 0);
end

mask_rcorr = false(height(T), 1);
if ismember('Win_len', T.Properties.VariableNames)
    idx_512 = (T.Win_len == 512);   mask_rcorr(idx_512) = T.Rcorr(idx_512) > 0.7;
    idx_1024 = (T.Win_len == 1024); mask_rcorr(idx_1024) = T.Rcorr(idx_1024) > 0.5;
    idx_2048 = (T.Win_len == 2048); mask_rcorr(idx_2048) = T.Rcorr(idx_2048) > 0.2;
    idx_4096 = (T.Win_len == 4096); mask_rcorr(idx_4096) = T.Rcorr(idx_4096) > 0.1;
    other = ~ismember(T.Win_len, [512, 1024, 2048, 4096]);
    mask_rcorr(other) = T.Rcorr(other) > Rcorr_Default;
else
    mask_rcorr = T.Rcorr > Rcorr_Default;
end

if ismember('Err_Az', T.Properties.VariableNames)
    mask_error = (T.Err_Az < Error_Threshold) & (T.Err_El < Error_Threshold);
else
    mask_error = true(height(T), 1);
end

final_idx = valid_geo & (~bad_region) & mask_rcorr & mask_error;
data_clean = T(final_idx, :);

%% ==================== 3. 提取数据与计算分形维数 ====================
idx_pos = data_clean.Start_loc >= Event1_Range(1) & data_clean.Start_loc <= Event1_Range(2);
data_pos = data_clean(idx_pos, :);

idx_neg = data_clean.Start_loc >= Event2_Range(1) & data_clean.Start_loc <= Event2_Range(2);
data_neg = data_clean(idx_neg, :);

% 计算 D 值 (注意：这里仍然使用全部数据计算，保证科学性，只在绘图时稀疏化)
if height(data_pos) > 50
    [D_pos, ~, ~] = calculateFractalDimension_BoxCount(data_pos.Azimuth, data_pos.Elevation);
else
    D_pos = NaN;
end
if height(data_neg) > 50
    [D_neg, ~, ~] = calculateFractalDimension_BoxCount(data_neg.Azimuth, data_neg.Elevation);
else
    D_neg = NaN;
end

%% ==================== 4. 绘图数据准备 ====================
if ~isempty(data_pos)
    t_pos_rel = (data_pos.Start_loc - min(data_pos.Start_loc)) / fs * 1e6;
end
if ~isempty(data_neg)
    t_neg_rel = (data_neg.Start_loc - min(data_neg.Start_loc)) / fs * 1e6;
end

%% ==================== 5. JGR 风格绘图 (指定Y轴范围并锁定比例) ====================
fig_width = 18; 
fig_height = 9;
figure('Units', 'centimeters', 'Position', [5, 5, fig_width, fig_height], 'Color', 'w');

% 布局参数
ax_bottom = 0.12; 
ax_height = 0.75; 
ax_width = 0.36;  
ax1_left = 0.07;  
ax2_left = 0.53;  
cb_left = 0.92;   

% 计算图框的长宽比 (Width/Height)
box_aspect_ratio = ax_width * fig_width / (ax_height * fig_height); 

font_name = 'Arial';
font_size = 10;
label_size = 11;
pt_size = 5; 

% --- (a) 左图：Pre-RS Negative Leader ---
ax1 = axes('Position', [ax1_left, ax_bottom, ax_width, ax_height]);
if ~isempty(data_neg)
    % 【关键修改】：定义绘图索引，从1开始，步长为2 (即隔一个取一个)
    plot_idx = 1:2:height(data_neg);
    
    % 使用 plot_idx 对 Azimuth, Elevation 和 Time 进行稀疏化采样
    scatter(data_neg.Azimuth(plot_idx), data_neg.Elevation(plot_idx), ...
            pt_size, t_neg_rel(plot_idx), 'filled', 'MarkerFaceAlpha', 0.8);
        
    text(0.05, 0.92, sprintf('D = %.2f', D_neg), 'Units', 'normalized', ...
        'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold', ...
        'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 3);
    
    % 设置指定的 YLim 并自动推算 XLim
    % 注意：这里建议仍使用 data_neg.Azimuth (全部数据) 来计算范围，以防稀疏化丢失边界点
    target_ylim = [0, 55]; 
    target_xlim = [165, 200]; 
    set_limits_fixed_ratio(ax1, target_ylim, target_xlim);
end
apply_jgr_style(ax1, font_name, font_size);
xlabel('Azimuth (°)', 'FontName', font_name, 'FontSize', label_size);
ylabel('Elevation (°)', 'FontName', font_name, 'FontSize', label_size);
title('(a) Pre-RS Negative Leader', 'FontName', font_name, 'FontSize', label_size, 'FontWeight', 'bold');


% --- (b) 右图：Post-RS Positive Leader ---
ax2 = axes('Position', [ax2_left, ax_bottom, ax_width, ax_height]);
if ~isempty(data_pos)
    scatter(data_pos.Azimuth, data_pos.Elevation, pt_size, t_pos_rel, 'filled', 'MarkerFaceAlpha', 0.8);
    text(0.05, 0.92, sprintf('D = %.2f', D_pos), 'Units', 'normalized', ...
        'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold', ...
        'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 3);
        
    % 设置指定的 YLim 并自动推算 XLim
    target_ylim = [40, 80]; 
    target_xlim = [150, 185]; 
    set_limits_fixed_ratio(ax2, target_ylim, target_xlim);
end
apply_jgr_style(ax2, font_name, font_size);
xlabel('Azimuth (°)', 'FontName', font_name, 'FontSize', label_size);
ylabel('Elevation (°)', 'FontName', font_name, 'FontSize', label_size); 
title('(b) Post-RS Positive Leader', 'FontName', font_name, 'FontSize', label_size, 'FontWeight', 'bold');


% --- Colorbar ---
cb = colorbar;
cb.Position = [cb_left, ax_bottom, 0.02, ax_height];
cb.Label.String = 'Relative Time (\mus)';
cb.Label.FontName = font_name;
cb.Label.FontSize = label_size;
colormap('jet'); 


%% ==================== 6. 辅助函数 ====================
function apply_jgr_style(ax, fname, fsize)
    set(ax, 'FontName', fname, 'FontSize', fsize, 'LineWidth', 1.0, ...
        'Box', 'on', 'TickDir', 'in', 'GridAlpha', 0.3, 'GridLineStyle', ':');
    grid on;
end

function set_limits_fixed_ratio(ax, target_ylim, target_xlim)
    ylim(ax, target_ylim);
    xlim(ax, target_xlim);
end

function [D, log_x, log_y] = calculateFractalDimension_BoxCount(az, el)
    if isempty(az), D=NaN; log_x=[]; log_y=[]; return; end
    az_rad = deg2rad(az); el_rad = deg2rad(el);
    x = cos(el_rad).*sin(az_rad); y = cos(el_rad).*cos(az_rad); z = sin(el_rad);
    points = [x, y, z];
    box_sizes = logspace(-1.8, -0.2, 8); 
    N = zeros(size(box_sizes));
    min_c = min(points,[],1); span = max(points,[],1)-min_c; span(span==0)=1;
    norm_p = (points - min_c)./span;
    for k=1:length(box_sizes)
        eps = box_sizes(k);
        idx = floor(norm_p * floor(1/eps)) + 1;
        N(k) = size(unique(idx,'rows'), 1);
    end
    valid = N>0;
    log_x = log(1./box_sizes(valid))'; log_y = log(N(valid))';
    if numel(log_x)<2, D=NaN; return; end
    p = polyfit(log_x, log_y, 1); D = p(1);
end