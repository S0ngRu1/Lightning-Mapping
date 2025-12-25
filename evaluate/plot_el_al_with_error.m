clear; clc; close all;

%% === 1. 参数设置 ===
% 文件名
filename = '..\2023\results\20230718175104_result_yld_3e8_6e8_window_1024_256_阈值4倍标准差_去零飘_20_80_hann_with_error.txt';

fs = 200e6;              % 采样率 200 MHz
Start_loc_Base = 5e8;    % 基准位置 (用于计算相对时间)
Error_Threshold = 0.5;   % 误差筛选阈值 (度)

%% === 2. 读取数据 ===
if ~isfile(filename)
    error('文件不存在: %s', filename);
end

% 使用 detectImportOptions 自动识别表头
opts = detectImportOptions(filename);
opts.VariableNamingRule = 'preserve'; % 保持原始列名
T = readtable(filename, opts);

% 检查是否存在误差列
if ~ismember('Err_Az', T.Properties.VariableNames) || ~ismember('Err_El', T.Properties.VariableNames)
    error('文件中未找到 Err_Az 或 Err_El 列，请检查文件格式。');
end

%% === 3. 数据筛选与处理 ===
% 1. 基本有效性筛选 (时间范围、仰角范围等)
valid_idx = T.Start_loc > Start_loc_Base & ...
            T.Start_loc < (Start_loc_Base + 8e7) & ...
            T.Elevation < 85 & ...
            T.Azimuth < 325 & ...
            T.Azimuth > 100 & ...
            T.Elevation > 25;

% 2. 核心筛选：误差小于阈值
valid_idx = valid_idx & (T.Err_Az < Error_Threshold) & (T.Err_El < Error_Threshold);

% 应用筛选
data = T(valid_idx, :);

if isempty(data)
    warning('没有符合筛选条件的数据点！');
    return;
end

% 3. 计算相对时间 (微秒 us)
data.Time_us = (data.Start_loc - Start_loc_Base) / fs * 1e6;

%% === 4. 绘图 ===

% 设置统一的绘图风格
figure('Units', 'pixels', 'Position', [100, 100, 1200, 900], 'Color', 'w');
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 颜色映射 (根据时间)
colors = data.Time_us; 

% --- 子图 1: 方位角误差随时间分布 ---
nexttile;
scatter(data.Time_us, data.Err_Az, 10, colors, 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
grid on; box on;
title('方位角误差随时间分布 (Azimuth Error vs Time)', 'FontWeight', 'bold');
ylabel('Azimuth Error (°)');
xlabel('Time (\mus)');
colorbar;
ylim([0, Error_Threshold * 1.1]); % 稍微留点余量

% --- 子图 2: 仰角误差随时间分布 ---
nexttile;
scatter(data.Time_us, data.Err_El, 10, colors, 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
grid on; box on;
title('仰角误差随时间分布 (Elevation Error vs Time)', 'FontWeight', 'bold');
ylabel('Elevation Error (°)');
xlabel('Time (\mus)');
colorbar;
ylim([0, Error_Threshold * 1.1]);

% --- 子图 3: Az-El 二维定位图 (带误差条) ---
% 这是一个合并的大图，占据下面两格
nexttile([1, 2]); 
hold on; box on; grid on;

% 1. 绘制误差条 (Error Bars)
% errorbar(x, y, y_neg, y_pos, x_neg, x_pos)
% 注意：为了避免图形太乱，误差条可以使用浅灰色，且不画所有点（如果点太多）
% 这里画所有点
e = errorbar(data.Azimuth, data.Elevation, ...
    data.Err_El, data.Err_El, ...  % Y轴误差 (上下)
    data.Err_Az, data.Err_Az, ...  % X轴误差 (左右)
    '.', 'Color', [0.6 0.6 0.6], 'CapSize', 0); % 灰色误差条，CapSize=0隐藏端点横线
e.Annotation.LegendInformation.IconDisplayStyle = 'off'; % 图例中不显示误差条

% 2. 绘制定位点 (带颜色)
scatter(data.Azimuth, data.Elevation, 3, colors, 'filled');

% 3. 装饰
title(sprintf('VHF辐射源二维定位'), 'FontWeight', 'bold');
xlabel('Azimuth (°)');
ylabel('Elevation (°)');
axis equal; % 保持角度比例一致
colormap('jet');
cb = colorbar;
ylabel(cb, 'Time (\mus)');

% 自动调整范围
margin = 2;
xlim([min(data.Azimuth)-margin, max(data.Azimuth)+margin]);
ylim([min(data.Elevation)-margin, max(data.Elevation)+margin]);

% 增加主标题
sgtitle('闪电定位误差分析与结果展示', 'FontSize', 14);