function [azimuth_deg, elevation_deg] = cart2sph_standard(vector)
%CART2SPH_STANDARD 将笛卡尔向量转换为标准方位角和仰角。
%   [azimuth_deg, elevation_deg] = CART2SPH_STANDARD(vector) 将
%   向量 [vx, vy, vz] (+X东, +Y北, +Z上) 转换为标准方位角（度，从+Y北向+X东）
%   和仰角（度，从XY平面向上）。

vx = vector(1);
vy = vector(2);
vz = vector(3);

% 计算水平距离
horizontal_distance = norm([vx, vy, 0]); % 或者 sqrt(vx^2 + vy^2)

% 计算仰角（从水平面向上）
% atan2d(y, x) 在 MATLAB 中是参数 (y, x)
elevation_deg = atan2d(vz, horizontal_distance);

% 计算方位角（从+Y向+X）
% atan2d(y, x)。我们需要从+Y向+X的角度，所以是 atan2d(vx, vy)
azimuth_deg = atan2d(vx, vy);

% 确保方位角在 0-360 度范围内
azimuth_deg = mod(azimuth_deg + 360, 360);

% 处理垂直向量和零向量的特殊情况
if horizontal_distance < 1e-10 % 使用一个小的容差
     if vz > 0
         elevation_deg = 90.0;
         azimuth_deg = 0.0; % 约定正上方的方位角为 0
     elseif vz < 0
         elevation_deg = -90.0;
         azimuth_deg = 0.0; % 约定正下方的方位角为 0
     else % 向量为零
         elevation_deg = NaN;
         azimuth_deg = NaN;
         disp('警告: cart2sph_standard: 输入零向量，方向未定义。');
     end
end

end