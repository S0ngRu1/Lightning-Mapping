clear; 
close all; 
clc;
%% ==================== 1. 参数设置 ====================
% --- 文件路径 ---
DATA_FILE = 'results\3d_win512_cost_cal_yld_chj_dtoa_3.6e8_5.6e8_with_initial_results.csv';
% --- 筛选参数 ---
START_LOC = 3.82e8; % 起始位置 (采样点)
END_LOC = 3.875e8; % 结束位置 (采样点)
UNIFIED_EPSILON = 120;     % 统一的 EPSILON (聚类半径)
% START_LOC = 3.66e8; % 起始位置 (采样点)
% END_LOC = 3.715e8; % 结束位置 (采样点)
% UNIFIED_EPSILON = 80;     % 统一的 EPSILON (聚类半径)
% --- 计算参数 ---
SAMPLING_RATE = 200e6;    % 200 MS/s
BIN_WIDTH_MS  = 1.5;      % 时间切片长度 (ms)
MIN_CLUSTER_POINTS = 2;   % DBSCAN 最小聚类点数
MIN_DT_STEPS = 1;         % 计算差分时至少跨越几个点
MAX_VELOCITY_THRESHOLD = 4e6; 
%% ==================== 2. 数据加载与基础筛选 ====================
fprintf('正在加载数据...\n');
all_match_results = readtable(DATA_FILE);
% 1. 基础物理条件筛选
conditions = ([all_match_results.dlta] < 20000) & ...
             ([all_match_results.x_dtoa] > -10000) & ([all_match_results.x_dtoa] < 6000) & ...
             ([all_match_results.y_dtoa] > -10000) & ([all_match_results.y_dtoa] < 0) & ...
             ([all_match_results.z_dtoa] > 0)      & ([all_match_results.z_dtoa] < 10000) & ...
             ([all_match_results.r_gccs] > 0.1) & ...
             (abs([all_match_results.R3_value]) < 10000);
% 2. 时间范围筛选
time_condition = ([all_match_results.yld_start_loc] >= START_LOC) & ...
                 ([all_match_results.yld_start_loc] <= END_LOC);
% 3. 空间区域剔除筛选
% 逻辑：找出在这个盒子里的点，取反（~）即为保留的点
in_exclusion_zone = ([all_match_results.x_dtoa] > -210) & ([all_match_results.x_dtoa] < 674) & ...
([all_match_results.y_dtoa] > -6082) & ([all_match_results.y_dtoa] < 4667) & ...
([all_match_results.z_dtoa] > 230) & ([all_match_results.z_dtoa] < 1200);
% 4. 综合筛选 (基础 & 时间 & 非排除区域)
filtered_indices = find(conditions & time_condition & ~in_exclusion_zone);
data = all_match_results(filtered_indices, :);
if isempty(data)
    error('筛选后无数据，请检查筛选条件。');
end
% 提取数据
raw_time_samples = [data.yld_start_loc]';
x_coords = [data.x_dtoa];
y_coords = [data.y_dtoa];
z_coords = [data.z_dtoa];
% 转换为相对时间 (ms)
start_time_sec = START_LOC / SAMPLING_RATE;
time_ms = (raw_time_samples / SAMPLING_RATE - start_time_sec) * 1000; 
fprintf('筛选后剩余点数: %d\n', length(time_ms));
%% ==================== 3. 按分段计算速度分布 ====================
fprintf('正在按 %.1f ms 分段计算速度...\n', BIN_WIDTH_MS);
% 定义时间网格 (Bins)
time_edges = floor(min(time_ms)) : BIN_WIDTH_MS : ceil(max(time_ms));
num_bins = length(time_edges) - 1;
% 存储用于绘图的数据
plot_time_groups = []; % X轴
plot_velocities  = []; % Y轴
plot_medians     = []; % 中位数
plot_median_time = []; % 中位数时间
for i = 1:num_bins
    t_start = time_edges(i);
    t_end   = time_edges(i+1);
    bin_center = (t_start + t_end) / 2;
    
    % 1. 提取当前区间内的所有点
    in_bin_idx = find(time_ms >= t_start & time_ms < t_end);
    
    if length(in_bin_idx) < MIN_CLUSTER_POINTS
        continue; % 点太少，跳过
    end
    
    % 提取坐标
    pts = [x_coords(in_bin_idx), y_coords(in_bin_idx), z_coords(in_bin_idx)];
    t_pts = time_ms(in_bin_idx);
    
    % 2. DBSCAN 去噪
    cluster_idx = dbscan(pts, UNIFIED_EPSILON, MIN_CLUSTER_POINTS);
    main_id = mode(cluster_idx(cluster_idx > 0));
    
    if isempty(main_id) || isnan(main_id)
        continue;
    end
    
    % 保留主类簇
    mask = (cluster_idx == main_id);
    pts_clean = pts(mask, :);
    t_clean   = t_pts(mask);
    
    if length(t_clean) < 2
        continue;
    end
    
    % 3. 排序
    [t_sorted, sort_order] = sort(t_clean);
    pts_sorted = pts_clean(sort_order, :);
    
    % 4. 计算该区间内的“瞬时速度”分布
    t_seconds = t_sorted / 1000; 
    
    % 计算位移
    dX = diff(pts_sorted(:,1));
    dY = diff(pts_sorted(:,2));
    dZ = diff(pts_sorted(:,3));
    dT = diff(t_seconds);
    
    % 注意：此处保留了您代码中的 dY^2 + dZ^2 (可能是为了计算YZ平面的投影速度)
    dist_seg = sqrt( dY.^2 + dZ.^2);
    
    % 简单的去噪 (去除时间重叠点)
    valid_seg = (dT > 1e-9); 
    v_segments = dist_seg(valid_seg) ./ dT(valid_seg);
    
    % 5. 存储数据
    if ~isempty(v_segments)
        v_segments = v_segments(v_segments < MAX_VELOCITY_THRESHOLD); 
        
        if ~isempty(v_segments)
            % 存入长数组用于 boxchart
            current_times = repmat(bin_center, length(v_segments), 1);
            
            plot_time_groups = [plot_time_groups; current_times];
            plot_velocities  = [plot_velocities; v_segments];
            
            % 计算中位数 (基于清洗后的数据)
            plot_medians = [plot_medians, median(v_segments)];
            plot_median_time = [plot_median_time, bin_center];
        end
    end
end
%% ==================== 4. 绘图与结果输出 ====================
fprintf('正在绘图...\n');
figure('Color', 'w', 'Position', [100, 100, 1000, 600]);

% 单位转换：10^6 m/s
y_data = plot_velocities / 1e6;   % 所有瞬时速度样本 (用于箱线图)
y_median = plot_medians / 1e6;    % 每个时间窗的中位数 (用于折线)

hold on;

if ~isempty(y_data)
    
    % === 【修改】同时计算两组统计数据 ===
    
    % 1. 基于“所有样本”的统计 (反映数据的离散程度)
    all_mean = mean(y_data);
    all_min  = min(y_data);
    all_max  = max(y_data);
    
    % 2. 基于“中位数趋势”的统计 (反映先导的主流速度)
    med_mean = mean(y_median);
    med_min  = min(y_median);
    med_max  = max(y_median);
    
    fprintf('==================================================\n');
    fprintf('【统计结果 A】基于所有瞬时样本 (All Samples):\n');
    fprintf('   平均值: %.2f x 10^6 m/s\n', all_mean);
    fprintf('   范围:   (%.2f ~ %.2f) x 10^6 m/s\n', all_min, all_max);
    fprintf('--------------------------------------------------\n');
    fprintf('【统计结果 B】基于中位数趋势 (Medians Only):\n');
    fprintf('   平均值: %.2f x 10^6 m/s\n', med_mean);
    fprintf('   范围:   (%.2f ~ %.2f) x 10^6 m/s\n', med_min, med_max);
    fprintf('   (建议论文描述使用此数据，受噪点影响小)\n');
    fprintf('==================================================\n');
    
    % =================================
    
    % 1. 绘制箱线图 (隐藏离群点)
    b = boxchart(plot_time_groups, y_data, ...
        'BoxFaceColor', 'white', ...
        'BoxEdgeColor', 'black', ...
        'WhiskerLineColor', 'black', ...
        'MarkerStyle', 'none', ...    
        'BoxWidth', BIN_WIDTH_MS * 0.8);
    
    % 2. 绘制中位数折线
    plot(plot_median_time, y_median, 'k-o', ...
        'LineWidth', 1.5, ...
        'MarkerFaceColor', 'w', ...
        'MarkerSize', 6);
else
    warning('没有符合条件(v < %.1e)的速度数据可绘图', MAX_VELOCITY_THRESHOLD);
end

hold off;

% --- 坐标轴设置 ---
xlabel('时间/ms', 'FontSize', 14);
ylabel('正先导通道速率/(10^6 m·s^{-1})', 'FontSize', 14);

% 安全设置 Y 轴范围 (基于中位数趋势微调，避免被极端值拉得太宽)
if ~isempty(y_median)
    ylim([min(y_median)*0.5, max(y_median)*1.5]);
end

% 强制 X 轴范围和刻度
if ~isempty(time_edges)
    xlim([min(time_edges), max(time_edges)]);
    % 自动调整刻度密度
    xticks(min(time_edges): (max(time_edges)-min(time_edges))/10 : max(time_edges)); 
end

set(gca, 'FontSize', 12, 'LineWidth', 1, 'TickDir', 'in');

if ~isempty(y_data)
    % 动态调整图例位置
    text_x = min(time_edges) + 0.02*(max(time_edges)-min(time_edges));
    y_lims = ylim; 
    text_y = y_lims(2) * 2; 
    text(text_x, text_y, '—○— 中位数', 'FontSize', 12);
end

grid off;