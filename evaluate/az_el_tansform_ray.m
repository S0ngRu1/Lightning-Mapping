%% ================================================================
%% PART 1: 几何解算 (Geometry Solver)
%% ================================================================
clc; clear; close all;

% === 1.1 站点参数 (Station Parameters) ===
Lat_A = 23.6392983; Lon_A = 113.5957504; Alt_A = 70.83; % Station A (VHF)
Lat_B = 23.621376;  Lon_B = 113.596166;  Alt_B = 46.28; % Station B (Optical)

% === 1.2 锚点观测值 ===
Az_A_meas = 179.624; El_A_meas = 48.1327;
Az_B_meas = 340.855; El_B_meas = 34.7503;

% 1. 计算 B 站相对位置 (ENU)
[xB, yB, zB] = geodetic2enu(Lat_B, Lon_B, Alt_B, Lat_A, Lon_A, Alt_A, wgs84Ellipsoid);
P_A = [0; 0; 0];      
P_B = [xB; yB; zB];   

% 2. 射线方向向量
[uA, vA, wA] = aer2enu_vec(Az_A_meas, El_A_meas); Dir_A = [uA; vA; wA];
[uB, vB, wB] = aer2enu_vec(Az_B_meas, El_B_meas); Dir_B = [uB; vB; wB];

% 3. 求解交会
LHS = [Dir_A, -Dir_B];
RHS = P_B - P_A;
sol = LHS \ RHS;
range_A_Anchor = abs(sol(1)); 
range_B_Anchor = abs(sol(2));

% 4. 计算 Bias
Point_On_A = P_A + range_A_Anchor * Dir_A;
Point_On_B = P_B + range_B_Anchor * Dir_B;
Target_XYZ = (Point_On_A + Point_On_B) / 2; 

[Az_A_ideal, ~, ~] = enu2aer_vec(Target_XYZ(1), Target_XYZ(2), Target_XYZ(3));
Bias_Rad = angdiff(deg2rad(Az_A_meas), deg2rad(Az_A_ideal));
Bias_Inst_Az = rad2deg(Bias_Rad);

% 5. 定义投影平面
Plane_Normal = [uA; vA; 0]; 
Plane_Normal = Plane_Normal / norm(Plane_Normal); 
Plane_D = dot(Target_XYZ, Plane_Normal);

%% ================================================================
%% PART 2: 数据变换 (Data Transformation)
%% ================================================================

% 读取数据
filename = '..\2023\results\20230718175104_result_yld_3e8_6e8_window_1024_256_阈值4倍标准差_去零飘_20_80_hann_with_error.txt';
Start_loc_Base = 3e8+9e7+9000000+1.2e6; 

result1 = readtable(filename);
logicalIndex = abs(result1.t123) < 1 & abs(result1.Rcorr) > 0.2 & ...
               result1.Start_loc < Start_loc_Base + 4.8e5 & ...
               result1.Start_loc > Start_loc_Base & ...
               result1.Elevation > 0;
data = result1(logicalIndex, :);

num_pts = height(data);
Az_B_new = zeros(num_pts, 1);
El_B_new = zeros(num_pts, 1);

for i = 1:num_pts
    % 1. 获取原始数据
    az_raw = data.Azimuth(i);
    el_raw = data.Elevation(i);
    
    % 2. 修正 Bias
    az_corr = mod(az_raw + Bias_Inst_Az, 360);
    
    % 3. 射线方向
    [u, v, w] = aer2enu_vec(az_corr, el_raw);
    Ray_Dir = [u; v; w]; 
    
    % 4. 射线-平面求交
    denominator = dot(Ray_Dir, Plane_Normal); 
    
    if abs(denominator) < 1e-4
        dist = range_A_Anchor;
    else
        dist = Plane_D / denominator;
    end
    
    if dist < 0, dist = range_A_Anchor; end
    
    % 5. 坐标变换
    P_Lightning = P_A + dist * Ray_Dir;
    Vec_B = P_Lightning - P_B;
    
    % 6. 反算 B 站视角 
    [az_b, el_b, ~] = enu2aer_vec(Vec_B(1), Vec_B(2), Vec_B(3));
    Az_B_new(i) = az_b;
    El_B_new(i) = el_b;
end

%% ================================================================
%% PART 3: JGR 格式绘图
%% ================================================================

fig_width_cm = 12;  
fig_height_cm = 9; 
figure('Units', 'centimeters', 'Position', [10, 10, fig_width_cm, fig_height_cm], ...
       'Color', 'w', 'PaperPositionMode', 'auto');

Start_loc = data.Start_loc;
c_values = (Start_loc - min(Start_loc)) / (max(Start_loc) - min(Start_loc));

% 1. 散点图
scatter(Az_B_new, El_B_new, 3, c_values, 'filled', 'MarkerFaceAlpha', 0.8);

% 2. 坐标轴格式化
set(gca, ...
    'FontName', 'Arial', ...
    'FontSize', 9, ...         
    'LineWidth', 1.0, ...      
    'Box', 'on', ...           
    'TickDir', 'out', ...      
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'XGrid', 'on', 'YGrid', 'on');

% 3. 标签
xlabel('Azimuth (deg)', 'FontName', 'Arial', 'FontSize', 10);
ylabel('Elevation (deg)', 'FontName', 'Arial', 'FontSize', 10);

% 4. 范围控制 (Axis Tight + Margin)
axis tight; 
x_lims = xlim; y_lims = ylim;
dx = diff(x_lims) * 2;
dy = diff(y_lims) * 0.05;
xlim([x_lims(1)-dx, x_lims(2)+dx]);
ylim([y_lims(1)-dy, y_lims(2)+dy]);

% 5. Colorbar
colormap('jet'); 
cb = colorbar;
cb.Label.String = 'Time (a.u.)'; 
cb.Label.FontName = 'Arial';
cb.Label.FontSize = 9;
cb.FontSize = 9;
cb.TickDirection = 'out';

fprintf('\n=== JGR 绘图完成 ===\n');
fprintf('Azimuth 范围已修改为 -180 到 180 度。\n');

%% === 辅助函数 ===
function [e, n, u] = aer2enu_vec(az, el)
    % 输入 Az 仍可以是 0-360，三角函数能自动处理
    el = deg2rad(el);
    az = deg2rad(az);
    u = sin(el);
    r = cos(el);
    e = r * sin(az);
    n = r * cos(az);
end

function [az, el, r] = enu2aer_vec(e, n, u)
    % 【关键修改】：计算结果改为 -180 到 180 度
    r = sqrt(e^2 + n^2 + u^2);
    el = rad2deg(asin(u/r));
    
    % atan2(Y,X) 本身返回 (-pi, pi]，对应 (-180, 180]
    % 原代码：az = mod(rad2deg(atan2(e, n)), 360);
    % 修改后：去掉 mod 即可
    az = rad2deg(atan2(e, n)); 
end

function d = angdiff(th1, th2)
    d = th2 - th1;
    d = mod(d + pi, 2*pi) - pi;
end