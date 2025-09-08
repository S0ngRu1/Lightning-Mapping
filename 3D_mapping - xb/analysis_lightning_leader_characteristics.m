clear; 
close all; 
clc;
%% ==================== 1. 参数设置与数据加载 ====================
% --- 用户需根据实际情况修改的参数 ---
DATA_FILE = 'result_yld_window512_th5n_all'; % 您的数据文件名
SAMPLING_RATE = 200e6;            % 数据采集卡采样率 (Hz), 200 MS/s
ASSUMED_HEIGHT = 1;            % 假设的放电平均高度 (米), 用于将角度转换成距离

% --- 速度计算参数 ---
TIME_WINDOW_VEL = 0.010;          % 速度计算的时间窗口 (秒), 10ms
% --- 列定义  ---
COL_START_LOC = 1;
COL_AZIMUTH   = 8;
COL_ELEVATION = 9;
COL_RCORR     = 10;
COL_T123      = 11;

% --- 加载数据 ---
fprintf('正在加载数据: %s\n', DATA_FILE);
try
    data = readmatrix(DATA_FILE);
catch
    error('数据文件加载失败，请检查文件名或文件格式是否正确。');
end


%% ==================== 2. 数据筛选 (新增部分) ====================
fprintf('正在进行数据筛选...\n');
% 根据用户提供的条件创建逻辑索引
logicalIndex = ...
    abs(data(:, COL_RCORR)) > 0.6 & ...      % 相关系数 > 0.6
    data(:, COL_START_LOC) < 4.1e8 & ...       % 采样点 < 6e8
    data(:, COL_START_LOC) > 3.8e8 & ...     % 采样点 > 3.8e8
    data(:, COL_ELEVATION) < 80 & ...        % 仰角 < 80度
    abs(data(:, COL_T123)) < 1;              % t123绝对值 < 1

% 应用筛选
data = data(logicalIndex, :);

fprintf('筛选后有效数据点数: %d\n', size(data, 1));

% 检查筛选后是否还有数据
if isempty(data)
    error('筛选后无数据点，请检查筛选条件或原始数据。');
end


% 提取所需列
time_samples = data(:, 1);  % 第1列: 起始采样点
azimuth_deg = data(:, 9);   % 第9列: 方位角 (度)
elevation_deg = data(:, 10);% 第10列: 仰角 (度)

%% ==================== 2. 数据预处理 ====================
% 将采样点转换为时间 (秒)，并以第一个点为时间零点
time_sec = (time_samples - time_samples(1)) / SAMPLING_RATE;

% 将角度(度)转换为弧度
azimuth_rad = deg2rad(azimuth_deg);
elevation_rad = deg2rad(elevation_deg);

% 将球面坐标(角度)投影到二维笛卡尔坐标系 (米)
% 采用切面投影法，假设放电发生在同一高度平面
% x: 东西方向, y: 南北方向, z: 垂直方向 (这里简化为二维)
% R = H / tan(el), x = R * cos(az), y = R * sin(az)
% 为简化，直接用 H*azimuth, H*elevation的近似也可，这里用更精确的投影
R_proj = ASSUMED_HEIGHT ./ tan(elevation_rad);
x_coords = R_proj .* cos(azimuth_rad);
y_coords = R_proj .* sin(azimuth_rad);

% 移除坐标中的非有限值 (可能由仰角接近0或90度引起)
valid_indices = isfinite(x_coords) & isfinite(y_coords);
time_sec = time_sec(valid_indices);
x_coords = x_coords(valid_indices);
y_coords = y_coords(valid_indices);

fprintf('数据预处理完成，共计 %d 个有效定位点。\n\n', length(time_sec));



%% ==================== 3. 调用函数进行计算与分析 ====================

% --- 计算发展速度 ---
fprintf('========== 开始计算发展速度 (时间窗口: %.0f ms) ==========\n', TIME_WINDOW_VEL * 1000);
[velocities, time_midpoints, avg_vel, std_vel] = calculate_velocity(time_sec, x_coords, y_coords, TIME_WINDOW_VEL);
fprintf('计算完成。\n');
fprintf('平均发展速度: %.2e ± %.2e m/s\n', avg_vel, std_vel);
plot(time_midpoints(2:9), velocities(2:9)/1e5, '-xblack', 'LineWidth', 1, ...
         'MarkerFaceColor', 'black', 'DisplayName', '负先导')

xlabel('时间 (ms)');
ylabel('发展速度 (\times10^5 m/s)');
hold on 
plot(time_midpoints(2:10), velocities(2:10)/1e5, ...
         'LineStyle', '-', ...
         'Marker', 'o', ...
         'Color', [0.5 0.5 0.5], ... % 使用RGB三元组定义中等灰色
         'LineWidth', 1, ...
         'DisplayName', '正先导')
legend('show', 'Location', 'best'); % 显示图例
grid on
function [velocities, time_midpoints, avg_vel, std_vel] = calculate_velocity(t, x, y, time_window)
% calculate_velocity: 计算先导在二维平面上的发展速度
% Inputs:
%   t            - 时间向量 (秒)
%   x, y         - 空间坐标向量 (米)
%   time_window  - 计算速度的时间窗口 (秒)
% Outputs:
%   velocities     - 每个时间窗口计算出的速度向量
%   time_midpoints - 每个速度对应的时间中点
%   avg_vel, std_vel - 平均速度和速度标准差

    % 初始化输出
    velocities = [];
    time_midpoints = [];

    % 计算总时长并确定分段数
    total_duration = t(end) - t(1);
    num_segments = floor(total_duration / time_window);

    if num_segments == 0
        warning('总时长小于一个时间窗口，无法进行分段速度计算。');
        avg_vel = NaN;
        std_vel = NaN;
        return;
    end

    % 循环遍历每个时间段
    for i = 1:num_segments
        t_start = t(1) + (i-1) * time_window;
        t_end = t_start + time_window;

        % 找到当前时间窗口内的所有数据点索引
        indices_in_window = find(t >= t_start & t < t_end);

        % 确保窗口内至少有2个点才能计算速度
        if length(indices_in_window) >= 2
            % 提取当前段的数据
            x_segment = x(indices_in_window);
            y_segment = y(indices_in_window);
            t_segment = t(indices_in_window);

            % 计算路径长度 (累加相邻点之间的欧氏距离)
            distances = sqrt(diff(x_segment).^2 + diff(y_segment).^2);
            path_length = sum(distances);
            
            % 计算时间跨度
            delta_t = t_segment(end) - t_segment(1);

            % 计算该段的速度
            if delta_t > 0
                current_velocity = path_length / delta_t;
                velocities = [velocities; current_velocity];
                time_midpoints = [time_midpoints; (i-2)*10];
            end
        end
    end

    % 计算平均值和标准差
    if ~isempty(velocities)
        avg_vel = mean(velocities);
        std_vel = std(velocities);
    else
        warning('所有时间段内的数据点均少于2个，无法计算任何速度。');
        avg_vel = NaN;
        std_vel = NaN;
    end
end

