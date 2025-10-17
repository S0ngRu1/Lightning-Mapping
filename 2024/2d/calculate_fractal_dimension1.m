%% ==================== 0. 初始化环境 ====================
clear; 
clc; 
close all;

%% ==================== 1. 参数设置与数据加载 ====================
% --- 用户需根据实际情况修改的参数 ---
DATA_FILE = 'result_yld_3.5e8_4e8_window_512_128_去零飘_滤波_加窗_阈值15_30_80.txt'; % 您的数据文件名

% --- 列定义  ---
COL_START_LOC = 1;
COL_AZIMUTH   = 8;
COL_ELEVATION = 9;
COL_RCORR     = 10;
COL_T123      = 11;

% --- 定义两个事件的时间窗口 (单位: 采样点) ---
event1_range = [3.65e8, 3.72e8]; % 例如，正先导
event2_range = [3.9e8, 4e8]; % 例如，负先导

% --- 加载原始数据 ---
fprintf('正在加载原始数据: %s\n', DATA_FILE);
try
    data_raw = readmatrix(DATA_FILE);
catch
    error('数据文件加载失败，请检查文件名或文件格式是否正确。');
end
fprintf('原始数据点数: %d\n', size(data_raw, 1));

%% ==================== 2. 分别处理和计算两个事件的分形维数 ====================

% --- 处理正先导 ---
fprintf('\n--- 正在处理正先导 (%.2e to %.2e) ---\n', event1_range(1), event1_range(2));
fractal_dim1 = NaN; log_x1 = []; log_y1 = []; % 初始化结果
try
    % 筛选正先导数据
    logicalIndex1 = ...
        abs(data_raw(:, COL_RCORR)) > 0.7 & ...
        data_raw(:, COL_START_LOC) > event1_range(1) & ...
        data_raw(:, COL_START_LOC) < event1_range(2) & ...
        data_raw(:, COL_ELEVATION) < 80 & ...
        abs(data_raw(:, COL_T123)) < 1;
    data1 = data_raw(logicalIndex1, :);
    fprintf('正先导筛选后有效数据点数: %d\n', size(data1, 1));
    
    if size(data1, 1) > 10 % 确保有足够的数据点
        az1 = data1(:, COL_AZIMUTH);
        el1 = data1(:, COL_ELEVATION);
        % 调用函数计算分形维数
        [fractal_dim1, log_x1, log_y1] = calculateFractalDimension_BoxCount(az1, el1);
        
    else
        fprintf('正先导数据点不足，跳过计算。\n');
    end
catch ME
    fprintf('处理正先导时出错: %s\n', ME.message);
end

% --- 处理负先导 ---
fprintf('\n--- 正在处理负先导 (%.2e to %.2e) ---\n', event2_range(1), event2_range(2));
fractal_dim2 = NaN; log_x2 = []; log_y2 = []; % 初始化结果
try
    % 筛选负先导数据
    logicalIndex2 = ...
        abs(data_raw(:, COL_RCORR)) > 0.6 & ...
        data_raw(:, COL_START_LOC) > event2_range(1) & ...
        data_raw(:, COL_START_LOC) < event2_range(2) & ...
        data_raw(:, COL_ELEVATION) < 80 & ...
        abs(data_raw(:, COL_T123)) < 1;
    data2 = data_raw(logicalIndex2, :);
    fprintf('负先导筛选后有效数据点数: %d\n', size(data2, 1));
    
    if size(data2, 1) > 10 % 确保有足够的数据点
        az2 = data2(:, COL_AZIMUTH);
        el2 = data2(:, COL_ELEVATION);
        % 调用函数计算分形维数
        [fractal_dim2, log_x2, log_y2] = calculateFractalDimension_BoxCount(az2, el2);
        
    else
        fprintf('负先导数据点不足，跳过计算。\n');
    end
catch ME
    fprintf('处理负先导时出错: %s\n', ME.message);
end

fprintf('正先导计算出的分形维数 D = %.3f\n', fractal_dim2);
fprintf('负先导计算出的分形维数 D = %.3f\n', fractal_dim1);
%% ==================== 3. 结果可视化 (精细布局版) ====================
fprintf('\n正在生成组合分析图...\n');
% --- 定义布局参数 ---
% 整体左右边距
margin_horiz = 0.1; 
% 顶部和底部图之间的垂直间距
gap_vert = 0.1;
% 顶部两个图之间的水平间距
gap_horiz = 0.08; 
% 底部图的高度
bottom_height = 0.35;
% 顶部图的高度
top_height = 0.38;

% --- 计算每个子图的位置 ---
% 下方子图的位置
bottom_width = 1 - 2*margin_horiz;
pos_bottom = [margin_horiz, 0.1, bottom_width, bottom_height];

% 上方子图的位置
top_width = (bottom_width - gap_horiz) / 2;
pos1 = [margin_horiz, pos_bottom(2) + bottom_height + gap_vert, top_width, top_height];
pos2 = [pos1(1) + top_width + gap_horiz, pos1(2), top_width, top_height];


% --- 左上子图: 正先导 空间分布 (方位角 vs 仰角) ---
axes('Position', pos1); % 使用手动计算的位置创建子图
if ~isempty(az1)
    plot(az1, el1, '.k', 'MarkerSize', 4);
    grid on;
    % 注意：对于角度图，axis equal 可能会过度拉伸，可手动设置范围
    title(sprintf('正先导 空间分布\n分形维数 D = %.3f', fractal_dim2), 'FontSize', 12);
    xlabel('方位角 (度)');
    ylabel('仰角 (度)');
    xlim([150 185]); % 保持您设定的范围
    ylim([30 78]);  % 保持您设定的范围
    set(gca, 'FontSize', 10);
else
    title('正先导: 无足够数据点');
end

% --- 右上子图: 负先导 空间分布 (方位角 vs 仰角) ---
axes('Position', pos2); % 使用手动计算的位置创建子图
if ~isempty(az2)
    plot(az2, el2, '.k', 'MarkerSize', 4);
    grid on;
    title(sprintf('负先导 空间分布\n分形维数 D = %.3f', fractal_dim1), 'FontSize', 12);
    xlabel('方位角 (度)');
    ylabel('仰角 (度)');
    xlim([130 200]); % 保持您设定的范围
    ylim([5 50]);   % 保持您设定的范围
    set(gca, 'FontSize', 10);
else
    title('负先导: 无足够数据点');
end

% --- 下方子图: 盒子计数法 log-log 对比图 ---
axes('Position', pos_bottom); % 使用手动计算的位置创建子图
hold on;
grid on;
box on;
if ~isnan(fractal_dim1), plot(log_x1, log_y1, 'bo', 'MarkerFaceColor', 'b'); p1 = polyfit(log_x1, log_y1, 1); fit1 = polyval(p1, log_x1); plot(log_x1, fit1, 'b-', 'LineWidth', 2); end
if ~isnan(fractal_dim2), plot(log_x2, log_y2, 'rs', 'MarkerFaceColor', 'r'); p2 = polyfit(log_x2, log_y2, 1); fit2 = polyval(p2, log_x2); plot(log_x2, fit2, 'r-', 'LineWidth', 2); end
title('盒子计数法对数-对数图', 'FontSize', 12);
xlabel('log(1/\epsilon)');
ylabel('log(N(\epsilon))');
legend_entries = {};
if ~isnan(fractal_dim1), legend_entries{end+1} = '正先导 数据点'; legend_entries{end+1} = sprintf('正先导 拟合 (D = %.3f)', fractal_dim2); end
if ~isnan(fractal_dim2), legend_entries{end+1} = '负先导 数据点'; legend_entries{end+1} = sprintf('负先导 拟合 (D = %.3f)', fractal_dim1); end
if ~isempty(legend_entries), legend(legend_entries, 'Location', 'northwest', 'FontSize', 10); end
set(gca, 'FontSize', 12);
hold off;

fprintf('绘图完成。\n');

%% ==================== 函数定义区 ====================
function [D, log_one_over_epsilon, log_N_epsilon] = calculateFractalDimension_BoxCount(azimuth_deg, elevation_deg)
    % calculateFractalDimension_BoxCount: 使用投影到单位球+三维盒子计数法计算分形维数
    % 输入: 
    %   azimuth_deg: 方位角向量 (角度制)
    %   elevation_deg: 仰角向量 (角度制)
    % 输出:
    %   D: 计算出的分形维数 (斜率)
    %   log_one_over_epsilon: 用于绘图的X轴数据 log(1/ε)
    %   log_N_epsilon: 用于绘图的Y轴数据 log(N(ε))
    
    % --- 1. 将角度坐标投影到单位球上 ---
    az_rad = deg2rad(azimuth_deg);
    el_rad = deg2rad(elevation_deg);
    
    x = cos(el_rad) .* sin(az_rad);
    y = cos(el_rad) .* cos(az_rad);
    z = sin(el_rad);
    
    points = [x, y, z];

    % --- 2. 三维盒子计数 ---
    % 定义一系列的盒子尺寸 ε (从0.01到1，共15个)
    box_sizes_epsilon = logspace(-2, 0, 5); 
    N_epsilon = zeros(size(box_sizes_epsilon));

    % 为了方便地将点分配到盒子里，我们先将所有点归一化到[0,1]的立方体内
    min_coords = min(points, [], 1);
    max_coords = max(points, [], 1);
    span = max_coords - min_coords;
    span(span == 0) = 1; % 防止分母为0
    
    normalized_points = (points - min_coords) ./ span;

    % 遍历每一种盒子尺寸
    for k = 1:length(box_sizes_epsilon)
        epsilon = box_sizes_epsilon(k);
        
        % 根据 ε 确定每个维度上的盒子数量
        num_boxes_dim = floor(1 / epsilon);
        if num_boxes_dim == 0, continue; end
        
        % 计算每个点落在哪一个盒子的索引中
        box_indices = floor(normalized_points * num_boxes_dim) + 1;
        
        % 找到所有被占据的盒子的唯一索引
        unique_boxes = unique(box_indices, 'rows');
        
        % 记录被占据的盒子数量
        N_epsilon(k) = size(unique_boxes, 1);
    end

    % --- 3. 准备 log-log 数据进行拟合 ---
    % 筛选出有效的数据点 (N(ε) > 0)
    valid_indices = N_epsilon > 0;
    
    log_one_over_epsilon = log(1 ./ box_sizes_epsilon(valid_indices))';
    log_N_epsilon = log(N_epsilon(valid_indices))';

    % --- 4. 线性拟合求解斜率 (即分形维数D) ---
    if numel(log_one_over_epsilon) < 2
        D = NaN;
        return;
    end
    
    p = polyfit(log_one_over_epsilon, log_N_epsilon, 1);
    D = p(1);
end