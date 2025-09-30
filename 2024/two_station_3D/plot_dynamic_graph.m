x_range = [-10000, 8000];
y_range = [-10000, 0];
z_range = [0, 10000];

conditions = ([all_match_results.chi_square_red] < 1500) & ...
    ([all_match_results.dlta] < 40000) & ...
    ([all_match_results.r_gccs] > 0.2) & ...
    (abs([all_match_results.R3_value]) < 1000 & ...
    ([all_match_results.yld_start_loc] >3.8e8) & ...
    ([all_match_results.yld_start_loc] < 5.5e8) );


% 一次性获取满足所有条件的索引
filtered_match_indices = find(conditions);

% 检查是否有数据通过筛选
if isempty(filtered_match_indices)
    error('没有数据点通过筛选条件，无法生成图像或视频。');
end

filtered_S_temp = all_S_results(filtered_match_indices, :);

% % 取优化前的原始三维结果
% initial_S_matrix = vertcat(all_match_results.S_initial_triangulation);
%
% filtered_S_temp = initial_S_matrix(filtered_match_indices, :);




filtered_match_result_temp = all_match_results(filtered_match_indices);
% 应用空间范围筛选
range_condition_s = filtered_S_temp(:,1) >= x_range(1) & filtered_S_temp(:,1) <= x_range(2) & ...
    filtered_S_temp(:,2) >= y_range(1) & filtered_S_temp(:,2) <= y_range(2) & ...
    filtered_S_temp(:,3) >= z_range(1) & filtered_S_temp(:,3) <= z_range(2);

filtered_S = filtered_S_temp(range_condition_s, :);
filtered_match_result = filtered_match_result_temp(range_condition_s);

% 提取坐标和时间数据
x = filtered_S(:, 1);
y = filtered_S(:, 2);
z = filtered_S(:, 3);
t = [filtered_match_result.yld_start_loc]'; % 明确这是时间向量

fprintf('数据筛选完成，共 %d 个点用于生成动画。\n', length(t));

[t_sorted, sort_idx] = sort(t);
x_sorted = x(sort_idx);
y_sorted = y(sort_idx);
z_sorted = z(sort_idx);

% 创建与时间演化对应的颜色图
num_points = length(t_sorted);
color_map = jet(num_points);

% --- 2.2 视频写入器设置 ---
video_filename = 'lightning_animation_short.mp4';
v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = 30;
v.Quality = 100;
open(v);

% --- 2.3 图形窗口和坐标轴初始化 ---
figure('Name', '动态闪电发展过程', 'Color', 'k', 'Position', [100, 100, 1000, 800]);
ax = gca;
ax.Color = 'k';
ax.GridColor = [0.7, 0.7, 0.7];
ax.GridAlpha = 0.5;

h_scatter = scatter3(nan, nan, nan, 2, nan, 'filled');
xlim(x_range);
ylim(y_range);
zlim(z_range);

xlabel('X (东 / m)', 'Color', 'w');
ylabel('Y (北 / m)', 'Color', 'w');
zlabel('Z (高度 / m)', 'Color', 'w');
title_handle = title('', 'Color', 'w', 'FontSize', 12);
ax.XColor = 'w'; ax.YColor = 'w'; ax.ZColor = 'w';
set(ax, 'DataAspectRatio', [1 1 1]);
daspect([1 1 1]);
view(130, 5);
grid on;
c = colorbar(ax);
c.Color = 'w';
ylabel(c, '归一化时间');
colormap(ax, color_map);
caxis([0, 1]);

fprintf('开始生成动画帧并写入视频...\n');

% --- 2.4 动画循环 ---
target_duration_sec = 15;
points_per_frame = ceil(num_points / (v.FrameRate * target_duration_sec));

points_per_frame = max(1, points_per_frame);


ts_ns = 5;

for frame_idx = 1:ceil(num_points / points_per_frame)
    end_point_idx = min(frame_idx * points_per_frame, num_points);
    current_indices = 1:end_point_idx;
    set(h_scatter, 'XData', x_sorted(current_indices), ...
        'YData', y_sorted(current_indices), ...
        'ZData', z_sorted(current_indices), ...
        'CData', color_map(current_indices, :));
    time_elapsed_ms = (t_sorted(end_point_idx) - t_sorted(1)) * ts_ns * 1e-6;
    set(title_handle, 'String', sprintf('闪电发展过程: 已用时间 %.2f ms (点: %d/%d)', ...
        time_elapsed_ms, end_point_idx, num_points));

    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% --- 2.5 清理工作 ---
close(v);
fprintf('视频已成功保存为: %s\n', video_filename);