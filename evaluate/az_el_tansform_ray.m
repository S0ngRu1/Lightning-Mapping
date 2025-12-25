clear; clc;

%% === 1. 站点坐标定义 (已填入您的数据) ===
% 您的 VHF 观测站 (Station A)
Site_A.Lat = 23.6392983; 
Site_A.Lon = 113.5957504; 
Site_A.Alt = 64.42;

% 光学观测站 (Station B)
Site_B.Lat = 23.621376; 
Site_B.Lon = 113.596166; 
Site_B.Alt = 56.930;

%% === 2. 模拟数据 (请替换为您真实的结果) ===
% 假设这是您算出的一组 VHF 结果
VHF_Az_List = [280; 285; 290]; % 方位角
VHF_El_List = [50;  55;  60];  % 仰角

%% === 3. 计算视差轨迹 ===
% 我们假设闪电可能发生的高度范围：2km 到 15km
Height_Range = 2000 : 100 : 15000; 

fprintf('正在计算从 VHF站 到 光学站 的视差转换...\n');

figure('Color','w', 'Position', [100, 100, 900, 600]);
hold on; grid on; axis equal;
xlabel('光学站方位角 (Azimuth @ Optical)');
ylabel('光学站仰角 (Elevation @ Optical)');
title('VHF 定位结果在光学视角的投影轨迹 (高度扫描法)');

% 颜色表，用于区分不同高度
cmap = jet(length(Height_Range));

% 遍历每一个 VHF 点
for i = 1:length(VHF_Az_List)
    curr_az_A = VHF_Az_List(i);
    curr_el_A = VHF_El_List(i);
    
    % 存储该点在 B 站视野中的轨迹
    Traj_Az_B = [];
    Traj_El_B = [];
    
    % --- 核心扫描循环 ---
    for h_idx = 1:length(Height_Range)
        h = Height_Range(h_idx);
        
        % 1. 在 A 站推算距离 (斜距 Range = 高差 / sin(El))
        delta_h = h - Site_A.Alt;
        if delta_h <= 0, continue; end % 忽略低于测站的点
        r_A = delta_h / sind(curr_el_A);
        
        % 2. 转为地心坐标 (ECEF)
        [x, y, z] = aer2ecef_manual(curr_az_A, curr_el_A, r_A, ...
                                    Site_A.Lat, Site_A.Lon, Site_A.Alt);
        
        % 3. 转为 B 站视角 (Az, El)
        [az_B, el_B, ~] = ecef2aer_manual(x, y, z, ...
                                          Site_B.Lat, Site_B.Lon, Site_B.Alt);
        
        Traj_Az_B(end+1) = az_B;
        Traj_El_B(end+1) = el_B;
    end
    
    % --- 绘图 ---
    % 画出这条线，代表这个 VHF 点在光学相机里可能出现的轨迹
    plot(Traj_Az_B, Traj_El_B, 'LineWidth', 2, 'DisplayName', sprintf('Pt %d (VHF: %.1f, %.1f)', i, curr_az_A, curr_el_A));
    
    % 标记几个关键高度点 (5km, 10km) 方便参照
    [~, idx_5km] = min(abs(Height_Range - 5000));
    [~, idx_10km] = min(abs(Height_Range - 10000));
    
    if idx_5km <= length(Traj_Az_B)
        plot(Traj_Az_B(idx_5km), Traj_El_B(idx_5km), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'g');
    end
    if idx_10km <= length(Traj_Az_B)
        plot(Traj_Az_B(idx_10km), Traj_El_B(idx_10km), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    end
end

legend show;
fprintf('计算完成。图中每条线代表一个 VHF 源点在光学视野中的可能位置。\n');
fprintf('绿色点 = 5km 高度，红色点 = 10km 高度。\n');


%% === 4. 底层坐标转换函数 (无需工具箱) ===

function [x, y, z] = aer2ecef_manual(az, el, r, lat0, lon0, h0)
    % 局部球坐标(Az,El,Range) -> 地心坐标(ECEF)
    az = deg2rad(az); el = deg2rad(el);
    % 1. 局部 ENU (East-North-Up)
    u = r .* sin(el);
    proj = r .* cos(el);
    n = proj .* cos(az);
    e = proj .* sin(az);
    
    % 2. ENU -> ECEF 旋转
    [x0, y0, z0] = lla2ecef_manual(lat0, lon0, h0);
    phi = deg2rad(lat0); lam = deg2rad(lon0);
    sl = sin(lam); cl = cos(lam); sp = sin(phi); cp = cos(phi);
    
    dx = -sl.*e - sp.*cl.*n + cp.*cl.*u;
    dy =  cl.*e - sp.*sl.*n + cp.*sl.*u;
    dz =          cp.*n     + sp.*u;
    
    x = x0 + dx; y = y0 + dy; z = z0 + dz;
end

function [az, el, r] = ecef2aer_manual(x, y, z, lat0, lon0, h0)
    % 地心坐标(ECEF) -> 局部球坐标(Az,El,Range)
    [x0, y0, z0] = lla2ecef_manual(lat0, lon0, h0);
    dx = x - x0; dy = y - y0; dz = z - z0;
    
    phi = deg2rad(lat0); lam = deg2rad(lon0);
    sl = sin(lam); cl = cos(lam); sp = sin(phi); cp = cos(phi);
    
    e = -sl.*dx + cl.*dy;
    n = -sp.*cl.*dx - sp.*sl.*dy + cp.*dz;
    u =  cp.*cl.*dx + cp.*sl.*dy + sp.*dz;
    
    r = sqrt(e.^2 + n.^2 + u.^2);
    el = rad2deg(asin(u ./ r));
    az = rad2deg(atan2(e, n));
    az = mod(az, 360);
end

function [x, y, z] = lla2ecef_manual(lat, lon, alt)
    % 经纬度转地心坐标 (WGS84)
    a = 6378137.0; f = 1 / 298.257223563; 
    b = a * (1 - f); e2 = 1 - (b^2 / a^2);
    lat = deg2rad(lat); lon = deg2rad(lon);
    N = a ./ sqrt(1 - e2 .* sin(lat).^2);
    
    x = (N + alt) .* cos(lat) .* cos(lon);
    y = (N + alt) .* cos(lat) .* sin(lon);
    z = (N * (1 - e2) + alt) .* sin(lat);
end