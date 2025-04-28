% --- 设置零点位置 ---
% 选择你要从哪个站点计算方位角和仰角
yld_sit = [0, 0, 0];       % 用户原始代码中的 yld_sit (行向量)
chj_sit = [1991, -7841.2, 0]; % 用户原始代码中的 chj_sit (行向量)

zero_point_position = yld_sit; % 选择 yld_sit 作为零点，你可以改成 chj_sit

% --- 使用 cart2sph_standard 计算方位角和仰角 ---

num_points = size(S_results, 1); % 定位结果的数量

% 预分配存储角度的数组
azimuths = zeros(num_points, 1);
elevations = zeros(num_points, 1);

% 遍历每一个定位结果点
for i = 1:num_points
    % 当前定位结果点 (确保是行向量)
    current_S = S_results(i, :);

    % 计算从零点位置指向当前定位结果点的向量
    % 确保输入到 cart2sph_standard 的是列向量
    direction_vector = current_S(:) - zero_point_position(:);

    % 调用 cart2sph_standard 函数计算方位角和仰角
    [azimuths(i), elevations(i)] = cart2sph_standard(direction_vector);
end

% --- 绘制二维图像 ---

% 方位角-仰角散点图
figure;
scatter(azimuths, elevations, 1, 'filled');
xlabel('方位角 (度)');
ylabel('仰角 (度)');
title(['从站点 ', num2str(zero_point_position), ' 看闪电的方位角 vs 仰角']);
grid on;
% 可以调整 X 轴范围以更好地展示方位角 (0-360度)
xlim([0, 360]);
% 可以调整 Y 轴范围以更好地展示仰角 (-90-90度)
ylim([0, 90]);