%% 新版数据结构
yld_sit = [0, 0, 0];
all_match_results = readtable('results\3d_win512_cost_cal_yld_chj_dtoa_3.6e8_5.6e8_with_initial_results_eval_dtoa.csv');
% 筛选条件
conditions = ([all_match_results.dlta] < 20000) & ...
             ([all_match_results.yld_start_loc] > 4.2e8) & ...
             ([all_match_results.yld_start_loc] < 5.6e8) & ...
             ([all_match_results.x_dtoa] > -10000) & ...
             ([all_match_results.x_dtoa] < 6000) & ...
             ([all_match_results.y_dtoa] > -10000) & ...
             ([all_match_results.y_dtoa] < 0) & ...
             ([all_match_results.z_dtoa] > 0) & ...
             ([all_match_results.z_dtoa] < 8000) & ...
             ([all_match_results.r_gccs] > 0.1) & ...
             (abs([all_match_results.R3_value]) < 10000);
filtered_match_indices = find(conditions);
filtered_match_result = all_match_results(filtered_match_indices, :);
% 处理颜色数据并归一化（增强对比度）
if ~isempty(filtered_match_result) &&  isnumeric([filtered_match_result.yld_start_loc])
    time_colors_raw = [filtered_match_result.yld_start_loc]';
    % 归一化并通过幂函数增强对比度（数值越大，对比度越强）
    time_colors = (time_colors_raw - min(time_colors_raw)) / (max(time_colors_raw) - min(time_colors_raw));
    time_colors = time_colors .^ 0.8;  % 降低幂次使颜色更鲜艳（0.5-0.8之间效果较好）
else
    disp('警告: filtered_match_result 为空，或 yld_start_loc 不可用/非数值类型。将按索引着色。');
    time_colors = (1:size(filtered_match_result,1))' / size(filtered_match_result,1);
    time_colors = time_colors .^ 0.8;  % 同样增强索引着色的对比度
end

marker_size = 5;  % 适当增大点大小，让颜色更显眼
x = [filtered_match_result.x_dtoa];
y = [filtered_match_result.y_dtoa]; 
z = [filtered_match_result.z_dtoa]; 

% --- 三维散点图 (鲜艳颜色风格 + 白色背景) ---
figure('Color', [1 1 1]); % 核心修改：设置白色背景

% 绘制三维散点图，提高点的不透明度让颜色更饱和
scatter3(x, y, z, marker_size, time_colors, 'filled', 'MarkerFaceAlpha', 0.8);

% 设置标题和轴标签（黑色字体，白色背景下更清晰）
xlabel('X (东)', 'FontSize', 12, 'Color', 'k'); 
ylabel('Y (北)', 'FontSize', 12, 'Color', 'k'); 
zlabel('Z (上)', 'FontSize', 12, 'Color', 'k'); 
title('辐射源三维定位空间分布', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');

% 设置坐标轴范围（保留原范围，确保数据显示正常）
% xlim([-10000, 6000]);
% ylim([-10000, 0]);
% zlim([0, 10000]);

% 坐标轴样式设置（黑色线条/文字，适配白色背景）
set(gca, ...
    'FontSize', 11, ...
    'LineWidth', 1.2, ...
    'Color', [1 1 1], ...  % 核心修改：坐标轴区域背景设为白色
    'XColor', 'k', ...     % 坐标轴刻度/线条设为黑色
    'YColor', 'k', ...
    'ZColor', 'k');

% 保留鲜艳颜色映射，确保数据颜色对比明显
colormap('jet');  % jet颜色映射高饱和，白色背景下显色清晰
% 备选高饱和度颜色映射：'hsv'、'hot'、'cool'、'rainbow'（可按需替换）
% colormap('hsv');  % 另一种鲜艳的选项

% 颜色条适配白色背景（黑色文字/线条）
h = colorbar;
ylabel(h, '归一化发展时间', 'FontSize', 11, 'Color', 'k');
set(h, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');  % 颜色条背景白、文字/刻度黑
caxis([0, 1]);  % 固定颜色范围以最大化对比度

% 网格和坐标轴比例（灰色网格避免干扰，适配白色背景）
grid on;
set(gca, ...
    'GridLineStyle', '--', ...
    'GridAlpha', 0.5, ...  % 适当提高网格透明度（白色背景下需更明显）
    'Box', 'on'); 
daspect([1 1 1]);  % 保持坐标轴比例一致，避免空间变形




%% 绘制三角定位数据
yld_sit = [0, 0, 0];
all_match_results = readtable('results\3d_win512_cost_cal_yld_chj_dtoa_3.6e8_5.6e8_with_initial_results.csv');
% 筛选条件
conditions = ([all_match_results.dlta] < 20000) & ...
             ([all_match_results.yld_start_loc] > 4.25e8) & ...
             ([all_match_results.yld_start_loc] < 5e8) & ...
             ([all_match_results.x_tri] > -10000) & ...
             ([all_match_results.x_tri] < 6000) & ...
             ([all_match_results.y_tri] > -10000) & ...
             ([all_match_results.y_tri] < 0) & ...
             ([all_match_results.z_tri] > 0) & ...
             ([all_match_results.z_tri] < 3500) & ...
             ([all_match_results.r_gccs] > 0.1) & ...
             (abs([all_match_results.R3_value]) < 10000);
filtered_match_indices = find(conditions);
filtered_match_result = all_match_results(filtered_match_indices, :);
% 处理颜色数据并归一化（增强对比度）
if ~isempty(filtered_match_result) &&  isnumeric([filtered_match_result.yld_start_loc])
    time_colors_raw = [filtered_match_result.yld_start_loc]';
    % 归一化并通过幂函数增强对比度（数值越大，对比度越强）
    time_colors = (time_colors_raw - min(time_colors_raw)) / (max(time_colors_raw) - min(time_colors_raw));
    time_colors = time_colors .^ 0.8;  % 降低幂次使颜色更鲜艳（0.5-0.8之间效果较好）
else
    disp('警告: filtered_match_result 为空，或 yld_start_loc 不可用/非数值类型。将按索引着色。');
    time_colors = (1:size(filtered_match_result,1))' / size(filtered_match_result,1);
    time_colors = time_colors .^ 0.8;  % 同样增强索引着色的对比度
end

marker_size = 20;  % 适当增大点大小，让颜色更显眼
x = [filtered_match_result.x_tri];
y = [filtered_match_result.y_tri]; 
z = [filtered_match_result.z_tri]; 

% --- 三维散点图 (鲜艳颜色风格 + 白色背景) ---
figure('Color', [1 1 1]); % 核心修改：设置白色背景

% 绘制三维散点图，提高点的不透明度让颜色更饱和
scatter3(x, y, z, marker_size, time_colors, 'filled', 'MarkerFaceAlpha', 0.8);

% 设置标题和轴标签（黑色字体，白色背景下更清晰）
xlabel('X (东)', 'FontSize', 12, 'Color', 'k'); 
ylabel('Y (北)', 'FontSize', 12, 'Color', 'k'); 
zlabel('Z (上)', 'FontSize', 12, 'Color', 'k'); 
title('辐射源三维定位空间分布', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');

% 设置坐标轴范围（保留原范围，确保数据显示正常）
% xlim([-10000, 6000]);
% ylim([-10000, 0]);
% zlim([0, 10000]);

% 坐标轴样式设置（黑色线条/文字，适配白色背景）
set(gca, ...
    'FontSize', 11, ...
    'LineWidth', 1.2, ...
    'Color', [1 1 1], ...  % 核心修改：坐标轴区域背景设为白色
    'XColor', 'k', ...     % 坐标轴刻度/线条设为黑色
    'YColor', 'k', ...
    'ZColor', 'k');

% 保留鲜艳颜色映射，确保数据颜色对比明显
colormap('jet');  % jet颜色映射高饱和，白色背景下显色清晰
% 备选高饱和度颜色映射：'hsv'、'hot'、'cool'、'rainbow'（可按需替换）
% colormap('hsv');  % 另一种鲜艳的选项

% 颜色条适配白色背景（黑色文字/线条）
h = colorbar;
ylabel(h, '归一化发展时间', 'FontSize', 11, 'Color', 'k');
set(h, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');  % 颜色条背景白、文字/刻度黑
caxis([0, 1]);  % 固定颜色范围以最大化对比度

% 网格和坐标轴比例（灰色网格避免干扰，适配白色背景）
grid on;
set(gca, ...
    'GridLineStyle', '--', ...
    'GridAlpha', 0.5, ...  % 适当提高网格透明度（白色背景下需更明显）
    'Box', 'on'); 
daspect([1 1 1]);  % 保持坐标轴比例一致，避免空间变形







% --- 4.2 二维投影图 (按时间着色) ---
figure; 
% XY 平面投影 (俯视图: 东-北)
scatter(x, y, marker_size, time_colors, 'filled');
xlabel('X (东)');
ylabel('Y (北)');
title('XY 平面投影 (东-北)');
grid on;
axis equal; 
colorbar;
colormap(gca, 'jet');

% XZ 平面投影 (侧视图: 东-上)
figure; 

scatter(x, z, marker_size, time_colors, 'filled'); 
xlabel('X (东)');
ylabel('Z (上)');
title('XZ 平面投影 (东-上)');
grid on;
axis equal;
colorbar;
colormap(gca, 'jet');

% YZ 平面投影 (前/后视图: 北-上)
figure; 

scatter(y, z, marker_size, time_colors, 'filled'); 
xlabel('Y (北)');
ylabel('Z (上)');
title('YZ 平面投影 (北-上)');
grid on;
axis equal;
colorbar;
colormap(gca, 'jet');

% --- 4.3极坐标系下 ---
theta = atan2(y, x); 
rho = sqrt(x.^2 + y.^2); 

figure;
polarscatter(theta, rho, marker_size, time_colors, 'filled');
title('极坐标系下');
colorbar;
colormap(gca, 'jet');


% --- 4.4 方位角 - 仰角 ---
num_points = size(filtered_match_result, 1);
azimuths = zeros(num_points, 1);
elevations = zeros(num_points, 1);

    for i = 1:num_points
    % 当前定位结果点 (确保是行向量)
    current_S = [filtered_match_result(i,:).x filtered_match_result(i,:).y filtered_match_result(i,:).z];
    % 计算从零点位置指向当前定位结果点的向量
    direction_vector = current_S(:) - yld_sit(:);
    % 调用 cart2sph_standard 函数计算方位角和仰角
    [azimuths(i), elevations(i)] = cart2sph_standard(direction_vector);
    end
    figure;
    scatter(azimuths, elevations, marker_size, time_colors, 'filled');
    xlabel('方位角 (度)');
ylabel('仰角 (度)');
    title(['从站点 ', num2str(yld_sit), ' 看闪电的方位角 vs 仰角']);
    grid on;
    xlim([0, 360]);
    xticks(0:40:360);
    ylim([0, 90]);
    yticks(0:10:90);
    colorbar;
    colormap(gca, 'jet');






%% 新版数据结构
yld_sit = [0, 0, 0];
all_match_results = readtable('results\3d_win512_cost_cal_yld_chj_dtoa_3.6e8_5.6e8_with_initial_results.csv');
% 筛选条件3d_win512_cost_cal_yld_chj_dtoa_3.6e8_5.6e8_with_initial_results_eval_dtoa
conditions = ([all_match_results.dlta] < 20000) & ...
             ([all_match_results.yld_start_loc] > 4.25e8) & ...
             ([all_match_results.yld_start_loc] < 5e8) & ...
             ([all_match_results.x_dtoa] > -10000) & ...
             ([all_match_results.x_dtoa] < 6000) & ...
             ([all_match_results.y_dtoa] > -10000) & ...
             ([all_match_results.y_dtoa] < 0) & ...
             ([all_match_results.z_dtoa] > 0) & ...
             ([all_match_results.z_dtoa] < 3500) & ...
             ([all_match_results.r_gccs] > 0.1) & ...
             (abs([all_match_results.R3_value]) < 10000);
filtered_match_indices = find(conditions);
filtered_match_result = all_match_results(filtered_match_indices, :);
% 处理颜色数据并归一化（增强对比度）
if ~isempty(filtered_match_result) &&  isnumeric([filtered_match_result.yld_start_loc])
    time_colors_raw = [filtered_match_result.yld_start_loc]';
    % 归一化并通过幂函数增强对比度（数值越大，对比度越强）
    time_colors = (time_colors_raw - min(time_colors_raw)) / (max(time_colors_raw) - min(time_colors_raw));
    time_colors = time_colors .^ 0.8;  % 降低幂次使颜色更鲜艳（0.5-0.8之间效果较好）
else
    disp('警告: filtered_match_result 为空，或 yld_start_loc 不可用/非数值类型。将按索引着色。');
    time_colors = (1:size(filtered_match_result,1))' / size(filtered_match_result,1);
    time_colors = time_colors .^ 0.8;  % 同样增强索引着色的对比度
end

marker_size = 20;  % 适当增大点大小，让颜色更显眼
x = [filtered_match_result.x_dtoa];
y = [filtered_match_result.y_dtoa]; 
z = [filtered_match_result.z_dtoa]; 

% --- 三维散点图 (鲜艳颜色风格 + 白色背景) ---
figure('Color', [1 1 1]); % 核心修改：设置白色背景

% 绘制三维散点图，提高点的不透明度让颜色更饱和
scatter3(x, y, z, marker_size, time_colors, 'filled', 'MarkerFaceAlpha', 0.8);

% 设置标题和轴标签（黑色字体，白色背景下更清晰）
xlabel('X (东)', 'FontSize', 12, 'Color', 'k'); 
ylabel('Y (北)', 'FontSize', 12, 'Color', 'k'); 
zlabel('Z (上)', 'FontSize', 12, 'Color', 'k'); 
title('辐射源三维定位空间分布', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');

% 设置坐标轴范围（保留原范围，确保数据显示正常）
% xlim([-10000, 6000]);
% ylim([-10000, 0]);
% zlim([0, 10000]);

% 坐标轴样式设置（黑色线条/文字，适配白色背景）
set(gca, ...
    'FontSize', 11, ...
    'LineWidth', 1.2, ...
    'Color', [1 1 1], ...  % 核心修改：坐标轴区域背景设为白色
    'XColor', 'k', ...     % 坐标轴刻度/线条设为黑色
    'YColor', 'k', ...
    'ZColor', 'k');

% 保留鲜艳颜色映射，确保数据颜色对比明显
colormap('jet');  % jet颜色映射高饱和，白色背景下显色清晰
% 备选高饱和度颜色映射：'hsv'、'hot'、'cool'、'rainbow'（可按需替换）
% colormap('hsv');  % 另一种鲜艳的选项

% 颜色条适配白色背景（黑色文字/线条）
h = colorbar;
ylabel(h, '归一化发展时间', 'FontSize', 11, 'Color', 'k');
set(h, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');  % 颜色条背景白、文字/刻度黑
caxis([0, 1]);  % 固定颜色范围以最大化对比度

% 网格和坐标轴比例（灰色网格避免干扰，适配白色背景）
grid on;
set(gca, ...
    'GridLineStyle', '--', ...
    'GridAlpha', 0.5, ...  % 适当提高网格透明度（白色背景下需更明显）
    'Box', 'on'); 
daspect([1 1 1]);  % 保持坐标轴比例一致，避免空间变形



%% 绘制三角定位数据
yld_sit = [0, 0, 0];
all_match_results = readtable('results\3d_win512_cost_cal_yld_chj_dtoa_3.6e8_5.6e8_with_initial_results_eval_dtoa.csv');
% 筛选条件
conditions = ([all_match_results.dlta] < 20000) & ...
             ([all_match_results.yld_start_loc] > 4.2e8) & ...
             ([all_match_results.yld_start_loc] < 5.6e8) & ...
             ([all_match_results.x_tri] > -10000) & ...
             ([all_match_results.x_tri] < 6000) & ...
             ([all_match_results.y_tri] > -10000) & ...
             ([all_match_results.y_tri] < 0) & ...
             ([all_match_results.z_tri] > 0) & ...
             ([all_match_results.z_tri] < 8000) & ...
             ([all_match_results.r_gccs] > 0.1) & ...
             (abs([all_match_results.R3_value]) < 10000);
filtered_match_indices = find(conditions);
filtered_match_result = all_match_results(filtered_match_indices, :);
% 处理颜色数据并归一化（增强对比度）
if ~isempty(filtered_match_result) &&  isnumeric([filtered_match_result.yld_start_loc])
    time_colors_raw = [filtered_match_result.yld_start_loc]';
    % 归一化并通过幂函数增强对比度（数值越大，对比度越强）
    time_colors = (time_colors_raw - min(time_colors_raw)) / (max(time_colors_raw) - min(time_colors_raw));
    time_colors = time_colors .^ 0.8;  % 降低幂次使颜色更鲜艳（0.5-0.8之间效果较好）
else
    disp('警告: filtered_match_result 为空，或 yld_start_loc 不可用/非数值类型。将按索引着色。');
    time_colors = (1:size(filtered_match_result,1))' / size(filtered_match_result,1);
    time_colors = time_colors .^ 0.8;  % 同样增强索引着色的对比度
end

x = [filtered_match_result.x_tri];
y = [filtered_match_result.y_tri]; 
z = [filtered_match_result.z_tri]; 

% --- 三维散点图 (鲜艳颜色风格 + 白色背景) ---
figure('Color', [1 1 1]); % 核心修改：设置白色背景

% 绘制三维散点图，提高点的不透明度让颜色更饱和
scatter3(x, y, z, marker_size, time_colors, 'filled', 'MarkerFaceAlpha', 0.8);

% 设置标题和轴标签（黑色字体，白色背景下更清晰）
xlabel('X (东)', 'FontSize', 12, 'Color', 'k'); 
ylabel('Y (北)', 'FontSize', 12, 'Color', 'k'); 
zlabel('Z (上)', 'FontSize', 12, 'Color', 'k'); 
title('辐射源三维定位空间分布', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');

% 设置坐标轴范围（保留原范围，确保数据显示正常）
% xlim([-10000, 6000]);
% ylim([-10000, 0]);
% zlim([0, 10000]);

% 坐标轴样式设置（黑色线条/文字，适配白色背景）
set(gca, ...
    'FontSize', 11, ...
    'LineWidth', 1.2, ...
    'Color', [1 1 1], ...  % 核心修改：坐标轴区域背景设为白色
    'XColor', 'k', ...     % 坐标轴刻度/线条设为黑色
    'YColor', 'k', ...
    'ZColor', 'k');

% 保留鲜艳颜色映射，确保数据颜色对比明显
colormap('jet');  % jet颜色映射高饱和，白色背景下显色清晰
% 备选高饱和度颜色映射：'hsv'、'hot'、'cool'、'rainbow'（可按需替换）
% colormap('hsv');  % 另一种鲜艳的选项

% 颜色条适配白色背景（黑色文字/线条）
h = colorbar;
ylabel(h, '归一化发展时间', 'FontSize', 11, 'Color', 'k');
set(h, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');  % 颜色条背景白、文字/刻度黑
caxis([0, 1]);  % 固定颜色范围以最大化对比度

% 网格和坐标轴比例（灰色网格避免干扰，适配白色背景）
grid on;
set(gca, ...
    'GridLineStyle', '--', ...
    'GridAlpha', 0.5, ...  % 适当提高网格透明度（白色背景下需更明显）
    'Box', 'on'); 
daspect([1 1 1]);  % 保持坐标轴比例一致，避免空间变形

