%%  2d静态图绘制
% --- 1. 数据准备  ---
filename = 'results\20240822165932_result_yld_3.65e8_5e8_window_256_64_阈值4倍标准差_去零飘_30_80_hann.txt';
Start_loc = 3.905e8;
% 2. 使用 readtable 函数读取数据
%    该函数会自动将第一行作为表头，并根据空格分隔各列
result1 = readtable(filename);
logicalIndex =  abs(result1.t123) < 0.5  & abs(result1.Rcorr) > 0.8 &  result1.Start_loc < Start_loc+2e6 & result1.Start_loc > Start_loc & result1.Elevation < 42 ;
filteredTable1 = result1(logicalIndex, :);

Start_loc = filteredTable1.Start_loc;
colorValues = (Start_loc - min(Start_loc)) / (max(Start_loc) - min(Start_loc)); % 归一化到 [0, 1]

% --- 2. 绘图  ---
% 设置figure的浅色背景
figure('Color', [1 1 1]); % figure背景设置为白色

% 使用 scatter 绘图，并应用尺寸和透明度优化
scatter(filteredTable1.Azimuth, filteredTable1.Elevation, ...
        10, ... % 尺寸
        colorValues, ...
        'filled', ...
        'MarkerFaceAlpha', 0.8); % 浅色背景下可适当提高透明度

% --- 3. 标签和标题优化 ---
% 设置标题和轴标签的颜色为深色
title('闪电VHF辐射源二维定位图', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
xlabel('方位角 (Azimuth / °)', 'FontSize', 12, 'Color', 'k');
ylabel('仰角 (Elevation / °)', 'FontSize', 12, 'Color', 'k');

% --- 4. 坐标轴和范围设置 ---
xlabel('方位角');
ylabel('仰角');
az_min = min(filteredTable1.Azimuth);
az_max = max(filteredTable1.Azimuth);
el_min = min(filteredTable1.Elevation);
el_max = max(filteredTable1.Elevation);
xlim([az_min, az_max]);  % 固定x轴范围（动态绘图不变化）
xticks('auto');          % 让Matlab在固定范围内自动生成规整刻度
ylim([el_min, el_max]);  % 固定y轴范围
yticks('auto');  

% 设置坐标轴的颜色和刻度字体颜色为深色
set(gca, ...
    'FontSize', 11, ...
    'LineWidth', 1.2, ...
    'Color', [1 1 1], ... % Axes背景色和figure背景色保持一致（白色）
    'XColor', [0.2 0.2 0.2], ... % X轴颜色（深灰色）
    'YColor', [0.2 0.2 0.2]);    % Y轴颜色（深灰色）

% --- 5. 颜色映射和颜色条优化 ---
colormap('jet'); % 保持专业的颜色映射
h = colorbar;

% 颜色条标签和刻度颜色为深色
ylabel(h, '归一化发展时间', 'FontSize', 11, 'Color', 'k');
set(h, 'Color', [0.2 0.2 0.2]); % 颜色条的刻度字体颜色（深灰色）

caxis([0, 1]); % 保持颜色范围

% --- 6. 网格和整体风格 ---
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3, 'Box', 'on'); % 浅色背景下降低网格透明度





%% 2d动态图绘制
filename = 'results\20240822165932_result_yld_3.6e8_3.8e8_window_512_128_阈值4倍标准差_去零飘_30_80_hann.txt';
Start_loc = 369650000;
% 2. 使用 readtable 函数读取数据
%    该函数会自动将第一行作为表头，并根据空格分隔各列
result1 = readtable(filename);
logicalIndex =  abs(result1.t123) < 0.5  & abs(result1.Rcorr) > 0.7 &  result1.Start_loc < Start_loc+1e5 & result1.Start_loc > Start_loc & result1.Elevation < 28 & result1.Azimuth < 184;
filteredTable1 = result1(logicalIndex, :);
Fs = 200e6; % 采样率 
target_viz_duration = 40; % (秒) 播放总时长
min_viz_pause = 5e-9; % (秒) 最小暂停时间

% 2. 计算真实时间和颜色
Start_loc = filteredTable1.Start_loc;
Start_loc_min = min(Start_loc);
Start_loc_max = max(Start_loc);
diff_start_loc = diff(Start_loc);
all_event_times_real = (Start_loc - Start_loc_min) / Fs;
all_colorValues = (Start_loc - Start_loc_min) / (Start_loc_max - Start_loc_min);
T_duration_real = max(all_event_times_real); % 闪电事件实际总时长

disp(['闪电事件实际持续时间: ', num2str(T_duration_real * 1e6), ' 微秒 (us)']);
disp(['将缩放到 ', num2str(target_viz_duration), ' 秒进行可视化回放...']);

% 3. 按时间顺序，对所有数据点进行排序
[sorted_times, sort_idx] = sort(all_event_times_real);
sorted_az = filteredTable1.Azimuth(sort_idx);
sorted_el = filteredTable1.Elevation(sort_idx);
sorted_colors = all_colorValues(sort_idx);

% 4. 计算点与点之间的真实时间间隔
% diff([0; V]) 会计算 t1-0, t2-t1, t3-t2 ...
time_gaps_real = diff([0; sorted_times]);

% 5. 计算缩放后的可视化暂停时间
% 缩放因子 = 目标总时长 / 真实总时长
scale_factor = target_viz_duration / T_duration_real;
viz_pauses = time_gaps_real * scale_factor;

% 6. 创建图形窗口
figure;
hold on;
grid on;
xlabel('方位角');
ylabel('仰角');
az_min = min(sorted_az);
az_max = max(sorted_az);
el_min = min(sorted_el);
el_max = max(sorted_el);
xlim([az_min, az_max]);  % 固定x轴范围（动态绘图不变化）
xticks('auto');          % 让Matlab在固定范围内自动生成规整刻度
ylim([el_min, el_max]);  % 固定y轴范围
yticks('auto');  

title('目标点动态呈现');
colormap('jet');
h_bar = colorbar;
ylabel(h_bar, '归一化起始位置');
caxis([0, 1]);

% 7. 逐点绘制，并按比例暂停
N = length(sorted_times);
disp(['开始回放... 总共 ', num2str(N), ' 个点。']);
tic; % 开始计时

for i = 1:N
    % (1) 获取当前点的暂停时间
    current_pause = viz_pauses(i);
    
    % (2) 执行暂停
    % 只有当暂停时间大于最小阈值时才真正暂停
    if current_pause > min_viz_pause
        pause(current_pause);
    end
    
    % (3) 绘制当前点
    scatter(sorted_az(i), sorted_el(i), 10, sorted_colors(i), 'filled');
    
    % (4) 刷新图像
    % 如果暂停时间太短，我们必须手动刷新，否则看不见
    if current_pause <= min_viz_pause
        drawnow; 
        % drawnow limitrate; % (如果点太多导致卡顿，可以尝试这个)
    end
end

viz_elapsed_time = toc; % 结束计时
hold off;
disp('回放完成。');
disp(['实际回放耗时: ', num2str(viz_elapsed_time), ' 秒 (理论目标: ', num2str(target_viz_duration), ' 秒)']);






%%  2d静态图绘制   大窗口长度
% --- 1. 数据准备  ---
filename = 'results\20240822165932_result_yld_3.6e8_5.6e8_window_1024_256_阈值4倍标准差_去零飘_30_80_hann.txt';
Start_loc = 370550000;
% 2. 使用 readtable 函数读取数据
%    该函数会自动将第一行作为表头，并根据空格分隔各列
result1 = readtable(filename);
logicalIndex =  abs(result1.t123) < 1  & abs(result1.Rcorr) > 0.6 &  result1.Start_loc < Start_loc+1e5 & result1.Start_loc > Start_loc & result1.Azimuth > 178 & result1.Azimuth < 184;
filteredTable1 = result1(logicalIndex, :);


Start_loc = filteredTable1.Start_loc;
colorValues = (Start_loc - min(Start_loc)) / (max(Start_loc) - min(Start_loc)); % 归一化到 [0, 1]

% --- 2. 绘图  ---
% 设置figure的浅色背景
figure('Color', [1 1 1]); % figure背景设置为白色

% 使用 scatter 绘图，并应用尺寸和透明度优化
scatter(filteredTable1.Azimuth, filteredTable1.Elevation, ...
        20, ... % 尺寸
        colorValues, ...
        'filled', ...
        'MarkerFaceAlpha', 0.8); % 浅色背景下可适当提高透明度

% --- 3. 标签和标题优化 ---
% 设置标题和轴标签的颜色为深色
title('闪电VHF辐射源二维定位图', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
xlabel('方位角 (Azimuth / °)', 'FontSize', 12, 'Color', 'k');
ylabel('仰角 (Elevation / °)', 'FontSize', 12, 'Color', 'k');

% --- 4. 坐标轴和范围设置 ---
xlabel('方位角');
ylabel('仰角');
az_min = min(filteredTable1.Azimuth);
az_max = max(filteredTable1.Azimuth);
el_min = min(filteredTable1.Elevation);
el_max = max(filteredTable1.Elevation);
xlim([az_min, az_max]);  % 固定x轴范围（动态绘图不变化）
xticks('auto');          % 让Matlab在固定范围内自动生成规整刻度
ylim([el_min, el_max]);  % 固定y轴范围
yticks('auto');  

% 设置坐标轴的颜色和刻度字体颜色为深色
set(gca, ...
    'FontSize', 11, ...
    'LineWidth', 1.2, ...
    'Color', [1 1 1], ... % Axes背景色和figure背景色保持一致（白色）
    'XColor', [0.2 0.2 0.2], ... % X轴颜色（深灰色）
    'YColor', [0.2 0.2 0.2]);    % Y轴颜色（深灰色）

% --- 5. 颜色映射和颜色条优化 ---
colormap('cool'); % 保持专业的颜色映射
h = colorbar;

% 颜色条标签和刻度颜色为深色
ylabel(h, '归一化发展时间', 'FontSize', 11, 'Color', 'k');
set(h, 'Color', [0.2 0.2 0.2]); % 颜色条的刻度字体颜色（深灰色）

caxis([0, 1]); % 保持颜色范围

% --- 6. 网格和整体风格 ---
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3, 'Box', 'on'); % 浅色背景下降低网格透明度


