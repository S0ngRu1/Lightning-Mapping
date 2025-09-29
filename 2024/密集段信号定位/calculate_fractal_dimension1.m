clear; 
close all; 
clc;
%% ==================== 1. 参数设置与数据加载 ====================
% --- 用户需根据实际情况修改的参数 ---
DATA_FILE = '20230718175104_result_yld_3e8_6e8_window_512_128_阈值4倍标准差_去零飘1_30_80_hann.txt';
SAMPLING_RATE = 200e6;            % 数据采集卡采样率 (Hz), 200 MS/s
ASSUMED_HEIGHT = 2500;            % 假设的放电平均高度 (米), 用于将角度转换成距离

% --- 列定义  ---
COL_START_LOC = 1;
COL_AZIMUTH   = 8;
COL_ELEVATION = 9;
COL_RCORR     = 10;
COL_T123      = 11;

% --- 定义两个事件的时间窗口 (单位: 采样点) ---
event1_range = [4.0e8, 4.2e8];
event2_range = [5.2e8, 5.4e8];

% --- 加载原始数据 ---
fprintf('正在加载原始数据: %s\n', DATA_FILE);
try
    data_raw = readmatrix(DATA_FILE);
catch
    error('数据文件加载失败，请检查文件名或文件格式是否正确。');
end
fprintf('原始数据点数: %d\n', size(data_raw, 1));


%% ==================== 2. 分别处理两个事件 ====================

% --- 处理事件1 ---
fprintf('\n--- 正在处理事件1 (%.2e to %.2e) ---\n', event1_range(1), event1_range(2));
fractal_dim1 = NaN; % 初始化结果
try
    % 筛选事件1数据
    logicalIndex1 = ...
        abs(data_raw(:, COL_RCORR)) > 0.6 & ...
        data_raw(:, COL_START_LOC) > event1_range(1) & ...
        data_raw(:, COL_START_LOC) < event1_range(2) & ...
        data_raw(:, COL_ELEVATION) < 80 & ...
        abs(data_raw(:, COL_T123)) < 1;
    data1 = data_raw(logicalIndex1, :);
    fprintf('事件1筛选后有效数据点数: %d\n', size(data1, 1));
    
    if size(data1, 1) > 10 % 确保有足够的数据点进行分析
        % 预处理事件1数据
        azimuth_deg1 = data1(:, COL_AZIMUTH);
        elevation_deg1 = data1(:, COL_ELEVATION);
        azimuth_rad1 = deg2rad(azimuth_deg1);
        elevation_rad1 = deg2rad(elevation_deg1);
        R_proj1 = ASSUMED_HEIGHT ./ tan(elevation_rad1);
        
        % 【修正】采用地理坐标系投影 (0度为正北, 顺时针)
        x_coords1 = R_proj1 .* sin(azimuth_rad1); % 东西方向 (East-West)
        y_coords1 = R_proj1 .* cos(azimuth_rad1); % 南北方向 (North-South)
        
        valid_indices1 = isfinite(x_coords1) & isfinite(y_coords1);
        x_coords1 = x_coords1(valid_indices1);
        y_coords1 = y_coords1(valid_indices1);
        
        % 计算事件1分形维数
        [fractal_dim1, ~, ~] = calculate_fractal_dimension(x_coords1, y_coords1);
        fprintf('事件1计算完成，分形维数 D = %.2f\n', fractal_dim1);
    else
        fprintf('事件1数据点不足，跳过计算。\n');
    end
end

% --- 处理事件2 ---
fprintf('\n--- 正在处理事件2 (%.2e to %.2e) ---\n', event2_range(1), event2_range(2));
fractal_dim2 = NaN; % 初始化结果
try
    % 筛选事件2数据
    logicalIndex2 = ...
        abs(data_raw(:, COL_RCORR)) > 0.6 & ...
        data_raw(:, COL_START_LOC) > event2_range(1) & ...
        data_raw(:, COL_START_LOC) < event2_range(2) & ...
        data_raw(:, COL_ELEVATION) < 80 & ...
        abs(data_raw(:, COL_T123)) < 1;
    data2 = data_raw(logicalIndex2, :);
    fprintf('事件2筛选后有效数据点数: %d\n', size(data2, 1));
    if size(data2, 1) > 10
        % 预处理事件2数据
        azimuth_deg2 = data2(:, COL_AZIMUTH);
        elevation_deg2 = data2(:, COL_ELEVATION);
        azimuth_rad2 = deg2rad(azimuth_deg2);
        elevation_rad2 = deg2rad(elevation_deg2);
        R_proj2 = ASSUMED_HEIGHT ./ tan(elevation_rad2);
        
        % 【修正】采用地理坐标系投影 (0度为正北, 顺时针)
        x_coords2 = R_proj2 .* sin(azimuth_rad2); % 东西方向 (East-West)
        y_coords2 = R_proj2 .* cos(azimuth_rad2); % 南北方向 (North-South)
        
        valid_indices2 = isfinite(x_coords2) & isfinite(y_coords2);
        x_coords2 = x_coords2(valid_indices2);
        y_coords2 = y_coords2(valid_indices2);

        % 计算事件2分形维数
        [fractal_dim2, ~, ~] = calculate_fractal_dimension(x_coords2, y_coords2);
        fprintf('事件2计算完成，分形维数 D = %.2f\n', fractal_dim2);
    else
        fprintf('事件2数据点不足，跳过计算。\n');
    end
end


%% ==================== 3. 结果可视化 (全新版本) ====================
fprintf('\n正在生成对比图...\n');
figure('Name', '双事件空间分布与分形维数', 'Position', [100, 100, 1200, 500]);

% --- 左子图: 事件1 空间分布图 ---
subplot(1, 2, 1);
if exist('x_coords1', 'var') && ~isempty(x_coords1)
    plot(x_coords1, y_coords1, '.k');
    grid on; axis equal;
    % 【修改】更新坐标轴标签和标题
    xlabel('东西方向 / m (East-West)');
    ylabel('南北方向 / m (North-South)');
    title(sprintf('负先导分形维数 (D = %.2f)', fractal_dim1));
    set(gca, 'FontSize', 12);
else
    title(sprintf('事件 1 (%.2e - %.2e)\n无足够数据点', event1_range(1), event1_range(2)));
    set(gca, 'FontSize', 12);
end

% --- 右子图: 事件2 空间分布图 ---
subplot(1, 2, 2);
if exist('x_coords2', 'var') && ~isempty(x_coords2)
    plot(x_coords2, y_coords2, '.k');
    grid on; axis equal;
    % 【修改】更新坐标轴标签和标题
    xlabel('东西方向 / m (East-West)');
    ylabel('南北方向 / m (North-South)');
    title(sprintf('正先导分形维数 (D = %.2f)', fractal_dim2));
    set(gca, 'FontSize', 12);
else
    title(sprintf('事件 2 (%.2e - %.2e)\n无足够数据点', event2_range(1), event2_range(2)));
    set(gca, 'FontSize', 12);
end

fprintf('绘图完成。\n');


%% ==================== 函数定义区 ====================
function [fractal_dim, box_counts, box_sizes] = calculate_fractal_dimension(x, y)
% calculate_fractal_dimension: 使用盒计数法计算二维点集的分形维数
    % --- 自适应设计BOX_SIZES参数 ---
    distances = sqrt(diff(x).^2 + diff(y).^2);
    avg_resolution = mean(distances(distances > 0));
    max_extent = max(max(x) - min(x), max(y) - min(y));
    r_min = ceil(avg_resolution); 
    if r_min == 0; r_min = 1; end
    r_max = floor(max_extent / 2);
    if r_max <= r_min; r_max = r_min + 100; end
    box_sizes = logspace(log10(r_min), log10(r_max), 20);
    
    % --- 盒计数法核心 ---
    box_counts = zeros(size(box_sizes));
    min_x = min(x); max_x = max(x);
    min_y = min(y); max_y = max(y);
    
    for i = 1:length(box_sizes)
        r = box_sizes(i);
        x_grid = min_x:r:max_x;
        y_grid = min_y:r:max_y;
        occupied_grid = false(length(y_grid), length(x_grid));
        ix = floor((x - min_x) / r) + 1;
        iy = floor((y - min_y) / r) + 1;
        linear_indices = sub2ind(size(occupied_grid), iy, ix);
        occupied_grid(linear_indices) = true;
        box_counts(i) = sum(occupied_grid(:));
    end
    
    valid_fit_indices = box_counts > 0;
    if sum(valid_fit_indices) < 2
        fractal_dim = NaN; return;
    end
    
    % --- 线性拟合计算维数 ---
    p = polyfit(log(box_sizes(valid_fit_indices)), log(box_counts(valid_fit_indices)), 1);
    fractal_dim = -p(1);
end