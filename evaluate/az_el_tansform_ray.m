clear; clc; close all;

%% ================================================================
%% PART 1: 几何诊断与参数计算 (Debug Geometry)
%% ================================================================

% === 1.1 站点参数 ===
Lat_A = 23.6392983; Lon_A = 113.5957504; Alt_A = 64.42; % 干涉仪(A)
Lat_B = 23.621376;  Lon_B = 113.596166;  Alt_B = 56.930; % 光学(B)

% === 1.2 锚点观测值 (请确保这是同一个闪电分支) ===
Az_A_meas = 202.651; El_A_meas = 50.3934;
Az_B_meas = 335.377; El_B_meas = 34.67;

fprintf('=== 步骤1: 几何诊断 ===\n');

% 1. 计算 B 站相对于 A 站的位置 (ENU坐标: 东, 北, 天)
[xB, yB, zB] = geodetic2enu(Lat_B, Lon_B, Alt_B, Lat_A, Lon_A, Alt_A, wgs84Ellipsoid);
P_A = [0; 0; 0];      
P_B = [xB; yB; zB];   

% 2. 计算两条射线的方向向量
[uA, vA, wA] = aer2enu_vec(Az_A_meas, El_A_meas); Dir_A = [uA; vA; wA];
[uB, vB, wB] = aer2enu_vec(Az_B_meas, El_B_meas); Dir_B = [uB; vB; wB];

% 3. 求解交会 (A + t*DirA = B + s*DirB)
LHS = [Dir_A, -Dir_B];
RHS = P_B - P_A;
sol = LHS \ RHS;
range_A = sol(1); % A 到目标的距离
range_B = sol(2); % B 到目标的距离

% 4. 几何合理性检查
fprintf('站点相对位置: B 在 A 的 [东: %.0fm, 北: %.0fm]\n', xB, yB);
fprintf('计算结果:\n');
fprintf('  A 站视线距离: %.2f km\n', range_A/1000);
fprintf('  B 站视线距离: %.2f km\n', range_B/1000);

if range_A < 0 || range_B < 0
    error('错误：计算出的距离为负数！这意味着两条射线方向发散，没有交点（或者交点在背后）。请检查输入的锚点角度是否正确匹配。');
end

% 5. 计算交会点与偏差
Point_On_A = P_A + range_A * Dir_A;
Point_On_B = P_B + range_B * Dir_B;
Target_XYZ = (Point_On_A + Point_On_B) / 2; % 空间最佳交点

% 计算系统偏差 (假设 B 是准的，强制把 A 校正到这个点)
[Az_A_ideal, ~, ~] = enu2aer_vec(Target_XYZ(1), Target_XYZ(2), Target_XYZ(3));
Bias_Inst_Az = Az_A_ideal - Az_A_meas;
fprintf('  >>> 干涉仪方位角修正量 (Bias): %+.2f°\n', Bias_Inst_Az);

% === 【诊断绘图】 画出两个站和射线，看看是否合理 ===
figure('Color','w', 'Position', [100,100,600,600]);
hold on; axis equal; grid on;
% 画站点
plot(0,0, 'b^', 'MarkerSize',10, 'LineWidth',2); text(0,0,' Station A');
plot(xB,yB, 'rs', 'MarkerSize',10, 'LineWidth',2); text(xB,yB,' Station B');
% 画射线 (画长一点，比如 10km)
L = 10000; 
quiver(0, 0, uA*L, vA*L, 'b--', 'MaxHeadSize',0, 'DisplayName', 'Ray A');
quiver(xB, yB, uB*L, vB*L, 'r--', 'MaxHeadSize',0, 'DisplayName', 'Ray B');
% 画交点
plot(Target_XYZ(1), Target_XYZ(2), 'kp', 'MarkerSize', 15, 'MarkerFaceColor','y');
text(Target_XYZ(1), Target_XYZ(2), ' Intersection');
xlabel('East (m)'); ylabel('North (m)'); 
title('俯视图：双站射线交会几何');
legend('Station A', 'Station B', 'Ray A', 'Ray B', 'Target');

%% ================================================================
%% PART 2: 数据变换 (Transform)
%% ================================================================
fprintf('\n=== 步骤2: 数据变换 ===\n');

% 读取数据
filename = '..\2023\results\20230718175104_result_yld_3e8_6e8_window_1024_256_阈值4倍标准差_去零飘_30_80_hann.txt';
Start_loc_Base = 3e8+9e7; 
if ~isfile(filename), error('文件不存在'); end
result1 = readtable(filename);

% 筛选
logicalIndex = abs(result1.t123) < 1 & abs(result1.Rcorr) > 0.2 & ...
               result1.Start_loc < Start_loc_Base + 3e7 & ...
               result1.Start_loc > Start_loc_Base & ...
               result1.Elevation > 0;
data = result1(logicalIndex, :);
fprintf('筛选后数据点: %d\n', height(data));

% === 核心变换循环 ===
num_pts = height(data);
Az_B_new = zeros(num_pts, 1);
El_B_new = zeros(num_pts, 1);

% 使用计算出的距离 range_A 进行投影
% 注意：这里假设所有点都位于 range_A 这个距离的球面上
for i = 1:num_pts
    % 1. 获取原始观测数据
    az_raw = data.Azimuth(i);
    el_raw = data.Elevation(i);
    
    % 2. 修正干涉仪的系统指北误差 (Bias)
    az_corr = az_raw + Bias_Inst_Az;
    
    % 3. 投影到 3D 空间 (相对于 A)
    % 公式: P = [e; n; u]
    [u, v, w] = aer2enu_vec(az_corr, el_raw);
    P_Lightning_A = [u; v; w] * range_A; 
    
    % 4. 坐标平移 (移到 B 站坐标系)
    % 向量 B->Lightning = (A->Lightning) - (A->B)
    Vec_B = P_Lightning_A - P_B;
    
    % 5. 反算 B 站视角
    [az_b, el_b, ~] = enu2aer_vec(Vec_B(1), Vec_B(2), Vec_B(3));
    Az_B_new(i) = az_b;
    El_B_new(i) = el_b;
end

%% ================================================================
%% PART 3: 绘图验证
%% ================================================================
figure('Color', [1 1 1], 'Position', [800, 100, 900, 700]);
Start_loc = data.Start_loc;
colorValues = (Start_loc - min(Start_loc)) / (max(Start_loc) - min(Start_loc)); 

scatter(Az_B_new, El_B_new, 15, colorValues, 'filled', 'MarkerFaceAlpha', 0.8);
hold on;

% 绘制锚点 (光学真值)
plot(Az_B_meas, El_B_meas, 'rp', 'MarkerSize', 20, 'LineWidth', 2, 'DisplayName', '光学锚点');
text(Az_B_meas, El_B_meas, ' Anchor(Optical)', 'FontSize', 12, 'Color', 'r');

title('转换结果 (光学 B 站视角)', 'FontSize', 16);
xlabel('Azimuth (°)'); ylabel('Elevation (°)');
grid on; colormap('jet'); colorbar;
legend('转换后的VHF点', '光学锚点真值');

% 自动范围
margin = 10;
xlim([min(Az_B_new)-margin, max(Az_B_new)+margin]);
ylim([min(El_B_new)-margin, max(El_B_new)+margin]);

%% === 辅助函数 (ENU坐标系: East, North, Up) ===
function [e, n, u] = aer2enu_vec(az, el)
    % 输入: 角度(度)
    % 输出: 单位向量 (东, 北, 天)
    % Az=0(北) -> e=0, n=1; Az=90(东) -> e=1, n=0
    el = deg2rad(el);
    az = deg2rad(az);
    
    u = sin(el);           % z (up)
    r = cos(el);           % 水平投影长度
    e = r * sin(az);       % x (east)
    n = r * cos(az);       % y (north)
end

function [az, el, r] = enu2aer_vec(e, n, u)
    % 输入: 向量 (东, 北, 天)
    % 输出: 角度(度)
    r = sqrt(e^2 + n^2 + u^2);
    el = rad2deg(asin(u/r));
    az = rad2deg(atan2(e, n)); % atan2(x, y) = atan2(East, North)
    az = mod(az, 360);
end