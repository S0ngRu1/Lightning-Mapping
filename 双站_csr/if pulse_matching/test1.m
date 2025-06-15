
num_points = size(filtered_S, 1); % 定位结果的数量

% 预分配存储角度的数组
azimuths = zeros(num_points, 1);
elevations = zeros(num_points, 1);
% 遍历每一个定位结果点
for i = 1:num_points
    % 当前定位结果点 (确保是行向量)
    current_S = filtered_S(i, :);

    % 计算从零点位置指向当前定位结果点的向量
    direction_vector = current_S(:) - yld_sit(:);

    % 调用 cart2sph_standard 函数计算方位角和仰角
    [azimuths(i), elevations(i)] = cart2sph_standard(direction_vector);
end
