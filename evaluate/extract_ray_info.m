%% === 1. 参数设置 (已根据您的相机参数修改) ===
% --- 相机硬件参数 ---
Img_Width = 1280;        % 图像宽度 (像素)
Img_Height = 1024;       % 图像高度 (像素)
Pixel_Size_mm = 4.8e-3;  % 像元尺寸: 4.8 um = 0.0048 mm
Focal_Length_mm = 5;     % 焦距: 5 mm

% --- 相机姿态参数 ---
% 仰角已知为 21度
Cam_Center_El = 21;      
% 方位角 (Azimuth) 仍需您根据实际拍摄方向确认
% 例如：如果相机朝正南拍，填180；朝正东拍，填90
Cam_Center_Az = 358.78;     

% --- 读取图像 ---
% 请修改为您实际的图片路径
img_filename = '0718_1751\LCI-HS_LRCHLS09_20230718-175103807.jpg'; 

if ~isfile(img_filename)
    warning('未找到图片文件，正在生成模拟图像用于演示...');
    [img_gray, ~, ~] = generate_mock_lightning(Img_Width, Img_Height);
else
    img_gray = imread(img_filename);
    if size(img_gray,3) > 1, img_gray = rgb2gray(img_gray); end
end

%% === 2. 图像处理：提取闪电通道 ===
% 2.1 阈值分割 (提取亮部)
% 简单的阈值分割，取最大亮度的 60%
threshold = 0.75 * double(max(img_gray(:))); 
binary_img = double(img_gray) > threshold;

% 去除噪点 (可选)
binary_img = bwareaopen(binary_img, 10);

% 2.2 骨架化 (Skeletonization)
skeleton_img = bwskel(binary_img);

% 2.3 获取像素坐标
[pixel_rows, pixel_cols] = find(skeleton_img);

% 坐标系转换：图像坐标(左上角为0,0) -> 图像物理坐标(中心为0,0)
% u: 水平方向 (列), v: 垂直方向 (行, 向上为正)
u = pixel_cols - Img_Width / 2;
v = (Img_Height / 2) - pixel_rows; 

%% === 3. 核心转换：像素 (u,v) -> 角度 (Az, El) ===
% 3.1 计算等效焦距 (以像素为单位)
% 公式：f_pixel = 物理焦距 / 像元尺寸
f_pixel = Focal_Length_mm / Pixel_Size_mm;

% (可选) 计算水平视场角 HFOV 用于验证
H_FOV_deg = 2 * atand((Img_Width * Pixel_Size_mm) / (2 * Focal_Length_mm));
fprintf('相机参数计算:\n - 等效焦距: %.2f pixels\n - 水平视场角: %.2f 度\n', f_pixel, H_FOV_deg);

% 3.2 局部相机坐标系向量 (Camera Frame)
% 向量 V_cam = [u, v, f]
num_pts = length(u);
V_cam = [u, v, ones(num_pts, 1) * f_pixel]'; % 3xN 矩阵

% 归一化向量
V_cam = V_cam ./ vecnorm(V_cam);

% 3.3 构建旋转矩阵 (Rotation Matrix)
% 定义相机视轴在世界系中的方向
az_rad = deg2rad(Cam_Center_Az);
el_rad = deg2rad(Cam_Center_El);

% 构建 ENU (东北天) 坐标系的旋转矩阵
% 假设相机无滚转角 (Roll=0)，且光轴对准 (Az, El)
% Forward (视轴)
F = [sin(az_rad)*cos(el_rad); cos(az_rad)*cos(el_rad); sin(el_rad)];
% Right (右轴) - 水平
R_vec = [sin(az_rad+pi/2); cos(az_rad+pi/2); 0]; 
% Up (上轴) - 正交于 F 和 R
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
opt_el = asin(zu ./ opt_range);       % 仰角 (弧度)
opt_az = atan2(xe, yn);               % 方位角 (弧度)

% 转为角度
Opt_El_Deg = rad2deg(opt_el);
Opt_Az_Deg = mod(rad2deg(opt_az), 360);

%% === 4. 绘图对比 ===
figure('Color', 'w', 'Position', [100, 100, 1200, 500]);

% 子图1：原始图像与提取结果
subplot(1, 2, 1);
imshow(img_gray); hold on;
[sk_r, sk_c] = find(skeleton_img);
plot(sk_c, sk_r, 'r.', 'MarkerSize', 1); 
title(['1. 提取的闪电通道 (' num2str(Img_Width) 'x' num2str(Img_Height) ')']);

% 子图2：转换后的 Az-El 坐标
subplot(1, 2, 2);
scatter(Opt_Az_Deg, Opt_El_Deg, 5, 'r', 'filled'); hold on;
grid on;
xlabel('Azimuth (°)'); ylabel('Elevation (°)');
title('2. 转换后的角度坐标 (Az-El)');
axis equal;

% 标记相机中心
plot(Cam_Center_Az, Cam_Center_El, 'k+', 'MarkerSize', 15, 'LineWidth', 2);
text(Cam_Center_Az, Cam_Center_El, ' Camera Center (21^\circ)');

fprintf('转换完成。\n');
fprintf('提取点数: %d\n', length(Opt_Az_Deg));
fprintf('方位角范围: %.1f ~ %.1f\n', min(Opt_Az_Deg), max(Opt_Az_Deg));
fprintf('仰角范围:   %.1f ~ %.1f\n', min(Opt_El_Deg), max(Opt_El_Deg));

%% === 辅助函数：生成模拟闪电图像 ===
function [img, az_t, el_t] = generate_mock_lightning(W, H)
    img = zeros(H, W);
    t = linspace(0, 1, 1000);
    % 模拟一条简单的闪电
    x_line = W/2 + 100 * sin(10*t); 
    y_line = H - (H-100) * t; 
    
    inds = sub2ind([H, W], round(y_line), round(x_line));
    img(inds) = 255;
    img = imgaussfilt(img, 2);
    img = uint8(mat2gray(img)*255);
    az_t = []; el_t = [];
end