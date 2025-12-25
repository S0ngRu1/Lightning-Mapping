clear; clc; close all;

%% === 1. 参数设置 (关键步骤：输入您的相机参数) ===
% 假设参数 (请根据您的高速相机实际情况修改)
Cam_H_FOV = 45;          % 水平视场角 (度)
Img_Width = 256;        % 图像宽度 (像素)
Img_Height = 512;       % 图像高度 (像素)

% 相机中心指向 (您当时相机对着哪里拍？)
% 如果不知道，通常需要找图像中的已知地标或恒星来反推
Cam_Center_Az = 180;     % 相机指向的中心方位角
Cam_Center_El = 45;      % 相机指向的中心仰角

% 读取图像 (这里生成模拟图像，实际使用请换成 imread('lightning.jpg'))
img_gray = imread('17_51_03\Img000014.jpg'); 
if size(img_gray,3)>1, img_gray = rgb2gray(img_gray); end
% [img_gray, true_Az, true_El] = generate_mock_lightning(Img_Width, Img_Height); % 生成模拟数据

%% === 2. 图像处理：提取闪电通道 ===
% 2.1 阈值分割 (提取亮部)
threshold = 0.6 * max(img_gray(:)); % 简单阈值，可换成 otus 或自适应
binary_img = img_gray > threshold;

% 2.2 骨架化 (Skeletonization)
% 将粗壮的闪电通道细化为单像素宽度的线
skeleton_img = bwskel(binary_img);

% 2.3 获取像素坐标
[pixel_rows, pixel_cols] = find(skeleton_img);

% 坐标系转换：图像坐标(左上角为0,0) -> 中心坐标(中心为0,0)
% u: 水平方向 (列), v: 垂直方向 (行, 向上为正)
u = pixel_cols - Img_Width / 2;
v = (Img_Height / 2) - pixel_rows; % 注意：图像y轴向下，物理仰角向上，需反转

%% === 3. 核心转换：像素 (u,v) -> 角度 (Az, El) ===
% 3.1 计算焦距 (以像素为单位)
% tan(FOV/2) = (Width/2) / f
f_pixel = (Img_Width / 2) / tand(Cam_H_FOV / 2);

% 3.2 局部相机坐标系向量 (Camera Frame)
% 假设相机坐标系：X右, Y上, Z前(光轴)
% 每个像素对应的向量 V_cam = [u, v, f]
num_pts = length(u);
V_cam = [u, v, ones(num_pts, 1) * f_pixel]'; % 3xN 矩阵

% 归一化向量
V_cam = V_cam ./ vecnorm(V_cam);

% 3.3 构建旋转矩阵 (Rotation Matrix)
% 将相机坐标系旋转到世界坐标系 (ENU: East-North-Up)
% 步骤：
% 1. 绕X轴旋转 -El_center (抬头)
% 2. 绕Z轴旋转 -Az_center (转头) -> 注意：方位角是顺时针，数学旋转是逆时针，需小心符号
% 这里采用标准导航坐标系转换：
% 定义相机视轴向量在世界系中的方向
az_rad = deg2rad(Cam_Center_Az);
el_rad = deg2rad(Cam_Center_El);

% 定义旋转矩阵 R (Camera -> World)
% 这是一个简化的旋转：先仰角 Pitch，再方位 Yaw
% R_el (绕X轴转)
R_el = [1, 0, 0;
        0, cos(el_rad), -sin(el_rad);
        0, sin(el_rad), cos(el_rad)];
    
% R_az (绕Z轴转，注意Az定义通常是北偏东，即Y偏X，这与标准数学定义不同)
% 为了简化，我们假设 Az=0 指向正北(Y)，Az=90 指向正东(X)
% 标准数学平面：X(东), Y(北). 
% 简单做法：先算出局部 Az/El 增量，再直接叠加 (小角度近似)
% 严谨做法：3D 旋转。

% === 严谨的 3D 旋转法 ===
% 1. 将局部 [u, v, f] 转为 [x, y, z] (Right, Up, Forward)
% 2. 旋转使其符合相机姿态
%    相机系: x(右), y(上), z(前)
%    世界系: E(东), N(北), U(天)

% 构建旋转矩阵：
% 相机z轴(视轴) 指向 (Az, El)
% 相机x轴(右) 指向 (Az+90, 0)
% 相机y轴(上) 指向 (Az, El+90)

% 视轴向量 (Forward)
F = [sin(az_rad)*cos(el_rad); cos(az_rad)*cos(el_rad); sin(el_rad)];
% 右轴向量 (Right) - 假设相机水平放置，没有Roll角
R_vec = [sin(az_rad+pi/2); cos(az_rad+pi/2); 0]; 
% 上轴向量 (Up) - 正交于 F 和 R
U_vec = cross(R_vec, F);

% 旋转矩阵 M = [R, U, F]
M_rot = [R_vec, U_vec, F]; 

% 3.4 转换到世界坐标系
V_world = M_rot * V_cam; % 3xN

% 3.5 从世界向量反算 Az, El
% V_world = [xe, yn, zu]
xe = V_world(1, :);
yn = V_world(2, :);
zu = V_world(3, :);

opt_range = sqrt(xe.^2 + yn.^2 + zu.^2);
opt_el = asin(zu ./ opt_range);       % 弧度
opt_az = atan2(xe, yn);               % atan2(x, y) = atan2(East, North) -> Azimuth

% 转为角度
Opt_El_Deg = rad2deg(opt_el);
Opt_Az_Deg = mod(rad2deg(opt_az), 360);

%% === 4. 绘图对比 ===
figure('Color', 'w', 'Position', [100, 100, 1200, 500]);

% 子图1：原始图像与提取结果
subplot(1, 2, 1);
imshow(img_gray); hold on;
[sk_r, sk_c] = find(skeleton_img);
plot(sk_c, sk_r, 'r.', 'MarkerSize', 1); % 绘制提取的骨架
title('1. 光学图像与提取的通道 (红色)');

% 子图2：转换后的 Az-El 坐标
subplot(1, 2, 2);
scatter(Opt_Az_Deg, Opt_El_Deg, 5, 'r', 'filled'); hold on;
grid on;
xlabel('Azimuth (°)'); ylabel('Elevation (°)');
title('2. 转换后的角度坐标 (Az-El)');
axis equal;
% 画出相机视场中心
plot(Cam_Center_Az, Cam_Center_El, 'k+', 'MarkerSize', 15, 'LineWidth', 2);
text(Cam_Center_Az, Cam_Center_El, ' Camera Center');

fprintf('转换完成。\n');
fprintf('提取点数: %d\n', length(Opt_Az_Deg));
fprintf('方位角范围: %.1f ~ %.1f\n', min(Opt_Az_Deg), max(Opt_Az_Deg));
fprintf('仰角范围:   %.1f ~ %.1f\n', min(Opt_El_Deg), max(Opt_El_Deg));


%% === 辅助函数：生成模拟闪电图像 ===
function [img, az_t, el_t] = generate_mock_lightning(W, H)
    img = zeros(H, W);
    
    % 模拟一条正弦波形状的闪电
    t = linspace(0, 1, 1000);
    x_line = 300 + 400 * t + 20 * sin(10*t); % X 轨迹
    y_line = 800 - 600 * t + 30 * cos(15*t); % Y 轨迹
    
    % 画到图像上
    inds = sub2ind([H, W], round(y_line), round(x_line));
    img(inds) = 255;
    % 高斯模糊模拟发光
    img = imgaussfilt(img, 2);
    img = uint8(mat2gray(img)*255);
    
    az_t = []; el_t = []; % 仅占位
end