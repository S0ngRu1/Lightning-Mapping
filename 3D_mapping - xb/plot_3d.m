x_range = [-10000, 6000];
y_range = [-10000, 0];
z_range = [0, 10000];

% 筛选条件
conditions = ([all_match_results.dlta] < 20000) & ...
             ([all_match_results.yld_start_loc] > 3.8e8) & ...
             ([all_match_results.yld_start_loc] < 3.9e8) & ...
             ([all_match_results.r_gccs] > 0.1) & ...
             (abs([all_match_results.R3_value]) < 10000);

% 获取满足条件的索引
filtered_match_indices = find(conditions);
filtered_S_temp = all_S_results(filtered_match_indices, :);
filtered_match_result_temp = all_match_results(filtered_match_indices); 

% 范围筛选
range_condition_s = filtered_S_temp(:,1) >= x_range(1) & filtered_S_temp(:,1) <= x_range(2) & ...
    filtered_S_temp(:,2) >= y_range(1) & filtered_S_temp(:,2) <= y_range(2) & ...
    filtered_S_temp(:,3) >= z_range(1) & filtered_S_temp(:,3) <= z_range(2);

filtered_S = filtered_S_temp(range_condition_s, :);
filtered_match_result = filtered_match_result_temp(range_condition_s);

% 处理颜色数据并归一化（增强对比度）
if ~isempty(filtered_match_result) && isfield(filtered_match_result, 'yld_start_loc') && isnumeric([filtered_match_result.yld_start_loc])
    time_colors_raw = [filtered_match_result.yld_start_loc]';
    % 归一化并通过幂函数增强对比度（数值越大，对比度越强）
    time_colors = (time_colors_raw - min(time_colors_raw)) / (max(time_colors_raw) - min(time_colors_raw));
    time_colors = time_colors .^ 0.8;  % 降低幂次使颜色更鲜艳（0.5-0.8之间效果较好）
else
    disp('警告: filtered_match_result 为空，或 yld_start_loc 不可用/非数值类型。将按索引着色。');
    time_colors = (1:size(filtered_S,1))' / size(filtered_S,1);
    time_colors = time_colors .^ 0.8;  % 同样增强索引着色的对比度
end

marker_size = 3;  % 适当增大点大小，让颜色更显眼
x = filtered_S(:, 1);
y = filtered_S(:, 2); 
z = filtered_S(:, 3); 

% --- 三维散点图 (鲜艳颜色风格) ---
figure('Color', [0.1 0.1 0.2]); % 深色背景更能凸显鲜艳颜色

% 绘制三维散点图，提高点的不透明度让颜色更饱和
scatter3(x, y, z, marker_size, time_colors, 'filled', 'MarkerFaceAlpha', 0.8);

% 设置标题和轴标签（白色字体增强对比）
xlabel('X (东)', 'FontSize', 12, 'Color', 'w'); 
ylabel('Y (北)', 'FontSize', 12, 'Color', 'w'); 
zlabel('Z (上)', 'FontSize', 12, 'Color', 'w'); 
title('辐射源三维定位空间分布', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'w');

% 设置坐标轴范围
xlim(x_range);
ylim(y_range);
zlim(z_range);

% 坐标轴样式设置
set(gca, ...
    'FontSize', 11, ...
    'LineWidth', 1.2, ...
    'Color', [0.1 0.1 0.2], ...
    'XColor', 'w', ...
    'YColor', 'w', ...
    'ZColor', 'w');

% 使用更鲜艳的颜色映射（可根据需要替换为其他高饱和度映射）
colormap('jet');  % jet颜色映射比parula更鲜艳
% 备选高饱和度颜色映射：'hsv'、'hot'、'cool'、'rainbow'
% colormap('hsv');  % 另一种鲜艳的选项

h = colorbar;
ylabel(h, '归一化发展时间', 'FontSize', 11, 'Color', 'w');
set(h, 'Color', 'w');
caxis([0, 1]);  % 固定颜色范围以最大化对比度

% 网格和坐标轴比例
grid on;
set(gca, ...
    'GridLineStyle', '--', ...
    'GridAlpha', 0.3, ...  % 降低网格透明度，避免干扰颜色
    'Box', 'on'); 
daspect([1 1 1]);


% 
% % --- 4.2 二维投影图 (按时间着色) ---
% figure; 
% % XY 平面投影 (俯视图: 东-北)
% scatter(x, y, marker_size, time_colors, 'filled');
% xlabel('X (东)');
% ylabel('Y (北)');
% title('XY 平面投影 (东-北)');
% grid on;
% axis equal; 
% colorbar;
% colormap(gca, 'cool');
% 
% % XZ 平面投影 (侧视图: 东-上)
% figure; 
% 
% scatter(x, z, marker_size, time_colors, 'filled'); 
% xlabel('X (东)');
% ylabel('Z (上)');
% title('XZ 平面投影 (东-上)');
% grid on;
% axis equal;
% colorbar;
% colormap(gca, 'cool');
% 
% % YZ 平面投影 (前/后视图: 北-上)
% figure; 
% 
% scatter(y, z, marker_size, time_colors, 'filled'); 
% xlabel('Y (北)');
% ylabel('Z (上)');
% title('YZ 平面投影 (北-上)');
% grid on;
% axis equal;
% colorbar;
% colormap(gca, 'cool');
% 
% % --- 4.3极坐标系下 ---
% theta = atan2(y, x); 
% rho = sqrt(x.^2 + y.^2); 
% 
% figure;
% polarscatter(theta, rho, marker_size, time_colors, 'filled');
% title('极坐标系下');
% colorbar;
% colormap(gca, 'cool');
% 
% 
% % --- 4.4 方位角 - 仰角 ---
% num_points = size(filtered_S, 1);
% azimuths = zeros(num_points, 1);
% elevations = zeros(num_points, 1);
% 
%     for i = 1:num_points
%     % 当前定位结果点 (确保是行向量)
%     current_S = filtered_S(i, :);
%     % 计算从零点位置指向当前定位结果点的向量
%     direction_vector = current_S(:) - yld_sit(:);
%     % 调用 cart2sph_standard 函数计算方位角和仰角
%     [azimuths(i), elevations(i)] = cart2sph_standard(direction_vector);
%     end
%     figure;
%     scatter(azimuths, elevations, marker_size, time_colors, 'filled');
%     xlabel('方位角 (度)');
% ylabel('仰角 (度)');
%     title(['从站点 ', num2str(yld_sit), ' 看闪电的方位角 vs 仰角']);
%     grid on;
%     xlim([0, 360]);
%     xticks(0:40:360);
%     ylim([0, 90]);
%     yticks(0:10:90);
%     colorbar;
%     colormap(gca, 'cool');