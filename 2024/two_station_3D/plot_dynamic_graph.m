all_match_results = readtable('results\3d_win512_cost_cal_yld_chj_dtoa3.6e8_4.0e8.csv');
% 筛选条件
conditions = ([all_match_results.dlta] < 20000) & ...
             ([all_match_results.yld_start_loc] > 3.6e8) & ...
             ([all_match_results.yld_start_loc] < 3.8e8) & ...
             ([all_match_results.x] > -10000) & ...
             ([all_match_results.x] < 6000) & ...
             ([all_match_results.y] > -10000) & ...
             ([all_match_results.y] < 0) & ...
             ([all_match_results.z] > 0) & ...
             ([all_match_results.z] < 10000) & ...
             ([all_match_results.r_gccs] > 0.1) & ...
             (abs([all_match_results.R3_value]) < 10000);
filtered_match_indices = find(conditions);
filtered_match_result = all_match_results(filtered_match_indices, :);
% 检查是否有数据通过筛选
if isempty(filtered_match_result)
    error('没有数据点通过筛选条件，无法生成图像或视频。');
end


% 提取坐标和时间数据
x = [filtered_match_result.x];
y = [filtered_match_result.y]; 
z = [filtered_match_result.z]; 
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
xlim([-10000 6000]);
ylim([-10000 0]);
zlim([0 10000]);

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