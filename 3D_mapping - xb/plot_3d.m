x_range = [-50000, 50000];
y_range = [-50000, 50000];
z_range = [0, 8000];
% 
conditions = ([all_match_results.chi_square_red] < 500) & ...
             ([all_match_results.yld_start_loc] > mapping_start_signal_loc) & ...
             ([all_match_results.yld_start_loc] < end_signal_loc) & ...
             ([all_match_results.r_gccs] > 0.1) & ...
             (abs([all_match_results.R3_value]) < 1000);


% conditions = ([all_match_results.dlta] < 20000) & ...
%              ([all_match_results.yld_start_loc] > mapping_start_signal_loc) & ...
%              ([all_match_results.yld_start_loc] < end_signal_loc) & ...
%              ([all_match_results.r_gccs] > 0.1) ;

% 一次性获取满足所有条件的索引
filtered_match_indices = find(conditions);
filtered_S_temp = all_S_results(filtered_match_indices, :);

% % 取优化前的原始三维结果
% initial_S_matrix = vertcat(all_match_results.S_initial_triangulation);
% 
% filtered_S_temp = initial_S_matrix(filtered_match_indices, :);

filtered_match_result_temp = all_match_results(filtered_match_indices); 
range_condition_s = filtered_S_temp(:,1) >= x_range(1) & filtered_S_temp(:,1) <= x_range(2) & ...
    filtered_S_temp(:,2) >= y_range(1) & filtered_S_temp(:,2) <= y_range(2) & ...
    filtered_S_temp(:,3) >= z_range(1) & filtered_S_temp(:,3) <= z_range(2);

filtered_S = filtered_S_temp(range_condition_s, :);
filtered_match_result = filtered_match_result_temp(range_condition_s);
if ~isempty(filtered_match_result) && isfield(filtered_match_result, 'yld_start_loc') && isnumeric([filtered_match_result.yld_start_loc])
    time_colors = [filtered_match_result.yld_start_loc]';
else
    disp('警告: filtered_match_result 为空，或 yld_start_loc 不可用/非数值类型。将按索引着色。');
    time_colors = (1:size(filtered_S,1))';
end

marker_size = 1;


x = filtered_S(:, 1);
y = filtered_S(:, 2); 
z = filtered_S(:, 3); 

% --- 4.1 三维散点图 (按时间着色) ---
figure;
%绘制三维散点图：x, y, z 为坐标，
scatter3(x, y, z, marker_size, time_colors, 'filled');
xlabel('X (东)'); 
ylabel('Y (北)'); 
zlabel('Z (上)'); 
title('源的三维空间分布');
grid on; 
axis equal;
colorbar; 
colormap(gca, 'cool');
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

% --- 4.4 方位角 - 仰角 ---
num_points = size(filtered_S, 1);
azimuths = zeros(num_points, 1);
elevations = zeros(num_points, 1);

    for i = 1:num_points
    % 当前定位结果点 (确保是行向量)
    current_S = filtered_S(i, :);
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
    colormap(gca, 'cool');