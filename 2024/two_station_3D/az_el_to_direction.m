function [x, y, z] = az_el_to_direction(azimuth, elevation)
    % 将方位角和仰角从度数转换为弧度
    azimuth_rad = deg2rad(azimuth);   % 方位角转换为弧度
    elevation_rad = deg2rad(elevation); % 仰角转换为弧度
    
    % 计算方向向量
    x = cos(elevation_rad) * cos(azimuth_rad);
    y = cos(elevation_rad) * sin(azimuth_rad);
    z = sin(elevation_rad);
end
