%% ========================================================================
%  后处理算法 ：筛选 -> 聚类 -> 可视化
% =========================================================================
clear; clc; close all;

%% 1. 加载原始定位结果
% -------------------------------------------------------------------------
filename = 'result_yld_3.8e8_4e8_window_512_128_去零飘_加窗_滤波_阈值15_30_80.txt';
if ~exist(filename, 'file')
    error('结果文件 "%s" 不存在，请检查文件名和路径。', filename);
end
fprintf('正在从 "%s" 读取原始定位结果...\n', filename);
results_table = readtable(filename);
fprintf('读取了 %d 个原始定位点。\n\n', height(results_table));

%% 2. 条件筛选
% -------------------------------------------------------------------------
logicalIndex = abs(results_table.t123) < 1 ...
             & abs(results_table.Rcorr) > 0.65 ...
             & results_table.Elevation < 80 ...
             & results_table.Start_loc < 4e8 ...
             & results_table.Start_loc > 3.8e8;

filtered_table = results_table(logicalIndex, :);
num_filtered_points = height(filtered_table);

fprintf('筛选完毕！剩下 %d / %d 个高质量定位点 (%.2f%%)。\n\n', ...
        num_filtered_points, height(results_table), 100 * num_filtered_points / height(results_table));

if num_filtered_points < 5
    fprintf('**警告：高质量定位点过少，无法进行有效的聚类聚合。**\n');
    return;
end

%% 3. 准备高质量数据用于聚类
% -------------------------------------------------------------------------
azimuths_deg = filtered_table.Azimuth;
elevations_deg = filtered_table.Elevation;

azimuths_rad = deg2rad(azimuths_deg);
elevations_rad = deg2rad(elevations_deg);
x = cos(elevations_rad) .* sin(azimuths_rad);
y = cos(elevations_rad) .* cos(azimuths_rad);
z = sin(elevations_rad);
location_data_3d = [x, y, z];


%% 4. 设置DBSCAN参数并执行聚类
% -------------------------------------------------------------------------
angular_distance_deg = 0.15;
epsilon = deg2rad(angular_distance_deg); 
MinPts = 1;

fprintf('正在对 %d 个高质量点进行DBSCAN聚类...\n', num_filtered_points);
[idx, ~] = dbscan(location_data_3d, epsilon, MinPts);
fprintf('聚类完成。\n\n');


%% 5. 处理并聚合聚类结果
% -------------------------------------------------------------------------
%% 5. 处理并聚合聚类结果 (新增：计算平均时间)
% -------------------------------------------------------------------------
fprintf('正在计算聚合后的中心点及其平均时间...\n');
cluster_indices = unique(idx(idx > 0));
num_clusters = length(cluster_indices);

if num_clusters == 0
    fprintf('**警告：未发现任何有效簇，请尝试放宽DBSCAN参数。**\n');
else
    fprintf('发现了 %d 个有效簇（聚合事件）。\n', num_clusters);
end

% 准备存储聚合后的位置和时间
aggregated_points_azel = zeros(num_clusters, 2);
aggregated_points_time = zeros(num_clusters, 1); % <--- 新增，用于存储平均时间

for i = 1:num_clusters
    current_cluster_id = cluster_indices(i);
    points_in_cluster_mask = (idx == current_cluster_id);
    
    % 计算中心位置 (矢量平均)
    vectors_in_cluster = location_data_3d(points_in_cluster_mask, :);
    mean_vector = mean(vectors_in_cluster, 1);
    mean_vector_normalized = mean_vector / norm(mean_vector);
    
    mx = mean_vector_normalized(1); my = mean_vector_normalized(2); mz = mean_vector_normalized(3);
    agg_az_rad = atan2(mx, my); agg_el_rad = asin(mz);
    
    agg_az_deg = rad2deg(agg_az_rad);
    if agg_az_deg < 0, agg_az_deg = agg_az_deg + 360; end
    agg_el_deg = rad2deg(agg_el_rad);
    
    aggregated_points_azel(i, :) = [agg_az_deg, agg_el_deg];

    % --- 【新增】计算这个簇的平均时间 ---
    mean_start_loc = mean(filtered_table.Start_loc(points_in_cluster_mask));
    aggregated_points_time(i) = mean_start_loc;
end


%% 6. 可视化一：详细的时空演化图
% -------------------------------------------------------------------------
fprintf('正在绘制图1：详细时空演化图...\n');
figure;
hold on;

% 计算原始点的颜色值
Start_loc_raw = filtered_table.Start_loc;
colorValues_raw = (Start_loc_raw - 3.8e8) / (4e8 - 3.8e8); 
color_range = [min(colorValues_raw), max(colorValues_raw)]; % 确定颜色范围

% 绘制所有筛选后的高质量点
scatter(filtered_table.Azimuth, filtered_table.Elevation, 1, colorValues_raw, 'filled');


% 设置图形属性
title('聚合前的放电事件', 'FontSize', 16);
xlabel('方位角 Azimuth (度)'); ylabel('俯仰角 Elevation (度)');
xlim([0, 360]); xticks(0:40:360); ylim([0, 90]); yticks(0:10:90); grid on;
colormap('cool');
h_cb = colorbar; ylabel(h_cb, '归一化时间');
if ~isempty(color_range), caxis(color_range); end % 统一颜色范围
hold off;


%% 7. 可视化二：仅含聚合结果的摘要图 (新增)
% -------------------------------------------------------------------------
if num_clusters > 0
    fprintf('正在绘制图2：仅含聚合事件的摘要图...\n');
    
    % --- 计算聚合点的颜色值 ---
    % 使用与第一张图完全相同的归一化方法，确保颜色意义一致
    colorValues_agg = (aggregated_points_time - 3.8e8) / (4e8 - 3.8e8);

    figure;
    hold on;

    % --- 绘制聚合后的点，颜色由其平均时间决定 ---
    scatter(aggregated_points_azel(:,1), aggregated_points_azel(:,2), ...
            1, colorValues_agg, 'filled'); 

    % --- 设置图形属性 ---
    title(sprintf('聚合后的 %d 个放电事件', num_clusters), 'FontSize', 16);
    xlabel('方位角 Azimuth (度)'); ylabel('俯仰角 Elevation (度)');
    xlim([0, 360]); xticks(0:40:360); ylim([0, 90]); yticks(0:10:90); grid on;
    
    % --- 设置与第一张图完全一致的颜色条 ---
    colormap('cool');
    h_cb2 = colorbar; ylabel(h_cb2, '归一化平均时间');
    if ~isempty(color_range), caxis(color_range); end % 使用和第一张图相同的颜色范围
    
    hold off;
end
fprintf('\n处理完毕。\n');