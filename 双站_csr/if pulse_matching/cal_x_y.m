% --- 主程序 ---
% CHJ 坐标
lat_chj = 23.5684654;
lon_chj = 113.6153785;
h_chj = -27; % 假设高度为 0

% YLD 坐标 (原点)
lat_yld = 23.6392983;
lon_yld = 113.5957504;
h_yld = 0; % 假设高度为 0

% 调用函数进行转换
[x_chj_rel_yld, y_chj_rel_yld, z_chj_rel_yld] = geodetic_to_local_enu(...
    lat_chj, lon_chj, h_chj, ...
    lat_yld, lon_yld, h_yld);

% 显示结果
fprintf('以 YLD 为原点 (0, 0, 0)，CHJ 的坐标为：\n');
fprintf('X (East) = %.4f 米\n', x_chj_rel_yld);
fprintf('Y (North)= %.4f 米\n', y_chj_rel_yld);
fprintf('Z (Up)   = %.4f 米 (假设两站高度均为0)\n', z_chj_rel_yld);

function [x_enu, y_enu, z_enu] = geodetic_to_local_enu(lat1, lon1, h1, lat0, lon0, h0)
% GEODETIC_TO_LOCAL_ENU 将经纬高转换为局部ENU坐标系
%   [x_enu, y_enu, z_enu] = GEODETIC_TO_LOCAL_ENU(lat1, lon1, h1, lat0, lon0, h0)
%   输入:
%       lat1, lon1, h1: 目标点的纬度(度)、经度(度)、高度(米)
%       lat0, lon0, h0: 原点的纬度(度)、经度(度)、高度(米)
%   输出:
%       x_enu, y_enu, z_enu: 目标点相对于原点的ENU坐标(米)

    % WGS84 椭球参数
    a = 6378137.0;         % 半长轴 (米)
    f = 1 / 298.257223563; % 扁率
    e_sq = 2*f - f^2;      % 第一偏心率的平方

    % 将经纬度从度转换为弧度
    phi1 = deg2rad(lat1);
    lambda1 = deg2rad(lon1);
    phi0 = deg2rad(lat0);
    lambda0 = deg2rad(lon0);

    % --- 将经纬高转换为 ECEF 坐标 ---
    function [X, Y, Z] = geodetic_to_ecef(phi, lambda, h)
        N = a / sqrt(1 - e_sq * sin(phi)^2); % 卯酉圈曲率半径
        X = (N + h) * cos(phi) * cos(lambda);
        Y = (N + h) * cos(phi) * sin(lambda);
        Z = (N * (1 - e_sq) + h) * sin(phi);
    end

    % 计算原点 (YLD) 的 ECEF 坐标
    [X0, Y0, Z0] = geodetic_to_ecef(phi0, lambda0, h0);

    % 计算目标点 (CHJ) 的 ECEF 坐标
    [X1, Y1, Z1] = geodetic_to_ecef(phi1, lambda1, h1);

    % 计算 ECEF 坐标差
    dX = X1 - X0;
    dY = Y1 - Y0;
    dZ = Z1 - Z0;

    % --- 将 ECEF 差值旋转到 ENU ---
    % 构建旋转矩阵 R
    R = [-sin(lambda0)           ,  cos(lambda0)           , 0;
         -sin(phi0)*cos(lambda0), -sin(phi0)*sin(lambda0), cos(phi0);
          cos(phi0)*cos(lambda0) ,  cos(phi0)*sin(lambda0) , sin(phi0)];

    % 应用旋转
    enu_coords = R * [dX; dY; dZ];

    x_enu = enu_coords(1); % East
    y_enu = enu_coords(2); % North
    z_enu = enu_coords(3); % Up
end

