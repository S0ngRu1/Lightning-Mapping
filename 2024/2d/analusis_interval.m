%% --- 1. 初始化和数据读取 ---
clear; clc; close all;

% 将您的数据保存到名为 'data.txt' 的文件中
filename = 'result_yld_3.8e8_4e8_window_512_128_去零飘_滤波_加窗_阈值15_30_80.txt';

% 使用 readtable 读取数据，它会自动处理表头和空格分隔符
try
    dataTable = readtable(filename);
catch ME
    error('无法读取文件 "%s"。请确保文件存在于当前目录，并且格式正确。', filename);
end

logicalIndex =  abs(dataTable.t123) < 1  & abs(dataTable.Rcorr) > 0.65 &  dataTable.Start_loc < 4e8 & dataTable.Start_loc > 3.8e8;
filteredTable1 = dataTable(logicalIndex, :);


%% --- 2. 定义时间区间并筛选数据 ---
interval1_start = 3.8e8;
interval1_end = 3.9e8;
interval2_start = 3.9e8;
interval2_end = 4.0e8;

% 使用逻辑索引筛选两个时间段内的数据
table_interval1 = filteredTable1(filteredTable1.Start_loc >= interval1_start & filteredTable1.Start_loc < interval1_end, :);
table_interval2 = filteredTable1(filteredTable1.Start_loc >= interval2_start & filteredTable1.Start_loc < interval2_end, :);

% 统计并打印数量
count1 = height(table_interval1);
count2 = height(table_interval2);

fprintf('--- 定位结果数量统计 ---\n');
fprintf('时间段 1 (%.1e - %.1e) 的定位点数量: %d\n', interval1_start, interval1_end, count1);
fprintf('时间段 2 (%.1e - %.1e) 的定位点数量: %d\n', interval2_start, interval2_end, count2);

%% --- 3. 分析负先导的发展特征 (时间段 2) ---
fprintf('\n--- 负先导发展特征分析 (%.1e - %.1e) ---\n', interval2_start, interval2_end);

if count2 > 1
    % MATLAB的sortrows可以对table进行排序
    table_interval2 = sortrows(table_interval2, 'Start_loc');

    % 使用 diff 函数计算相邻定位点之间的时间差
    time_diffs = diff(table_interval2.Start_loc);

    figure; % 创建一个新的图形窗口
    h = histogram(time_diffs);

    % 为了更好地观察分布，特别是当数值范围很大时，使用对数坐标轴
    set(gca, 'XScale', 'log');

    % 添加标题和标签，使图形更清晰
    title('Distribution of Time Differences to Determine Gap Threshold');
    xlabel('Time Difference (sampling points) - Log Scale');
    ylim([0,300])
    ylabel('Frequency (Count)');
    grid on; % 添加网格线

    % 定义阈值来区分“阶梯内部发展”和“阶梯间断”
    % 如果时间差大于这个值，我们认为是一次间断。这个值可以根据数据特征调整。
    gap_threshold = 39000; % 单位：采样点

    % 找出所有的间断间隔
    gaps = time_diffs(time_diffs > gap_threshold);
    % 找到间断发生的位置索引
    gap_indices = find(time_diffs > gap_threshold);

    % 阶梯的起始点索引：第一个点，以及每个间断后的第一个点
    step_start_indices = [1; gap_indices + 1];

    % 阶梯的结束点索引：每个间断前的点，以及最后一个点
    step_end_indices = [gap_indices; height(table_interval2)];


    fprintf('\n--- 各发展阶段的起始 Start_loc ---\n');

    % 直接使用 step_start_indices 从 table_interval2 中提取每个阶段的第一个 Start_loc 值
    % 这是一个矢量化操作，无需循环
    each_step_start_loc = table_interval2.Start_loc(step_start_indices);

    % 在命令行窗口中显示结果
    fprintf('每个阶段的起始位置 (Start_loc) 如下:\n');
    disp(each_step_start_loc);
end



%% --- 按真实时间发展的动态图绘制 ---
% 1. 为每个定位点分配一个阶段ID
num_steps = numel(step_start_indices);
table_interval2.StepID = zeros(height(table_interval2), 1); 
for i = 1:num_steps
    start_idx = step_start_indices(i);
    end_idx = step_end_indices(i);
    table_interval2.StepID(start_idx:end_idx) = i;
end

% 2. 为每个阶段ID创建一种独特的颜色
step_colors = lines(num_steps);

% --- 参数定义 ---
fs = 200e6;         % 采样频率 (Hz)
speed_factor = 500;  % 动画速度的放大倍数
plot_batch_size = 20; % 定义每次绘制的点数（批次大小）

% --- 图形初始化 ---
figure;
hold on;
grid on;
box on; 
xlabel('方位角 (Azimuth)');
xlim([120, 200]);
xticks(120:20:200);
ylabel('仰角 (Elevation)');
ylim([10, 70]);
yticks(10:10:70);
title('闪电发展阶段动态呈现 ');

% 创建图例
legend_handles = gobjects(num_steps, 1);
legend_entries = cell(num_steps, 1);
for i = 1:num_steps
    legend_handles(i) = scatter(NaN, NaN, 50, step_colors(i,:), 'filled');
    legend_entries{i} = sprintf('阶段 %d', i);
end
legend(legend_handles, legend_entries, 'Location', 'northwest');


for i = 1:plot_batch_size:height(table_interval2)
    % 1. 确定当前批次的起止行号
    start_idx = i;
    end_idx = min(i + plot_batch_size - 1, height(table_interval2));
    
    % 2. 提取当前批次的数据
    current_batch = table_interval2(start_idx:end_idx, :);
    
    % 如果批次为空，则跳过
    if isempty(current_batch)
        continue;
    end
    
    % 3. 绘制当前批次的所有点
    batch_step_ids = current_batch.StepID;
    valid_ids = batch_step_ids > 0;
    if any(valid_ids)
        batch_colors = step_colors(batch_step_ids(valid_ids), :);
        scatter(current_batch.Azimuth(valid_ids), current_batch.Elevation(valid_ids), 5, batch_colors, 'filled', 'HandleVisibility', 'off');
    end
    
    % 4. 刷新画布
    drawnow;
    
    % 5. 计算这一批次所跨越的真实时间，并据此计算暂停时长
    % 获取批次中第一个点和最后一个点的Start_loc
    batch_start_time = current_batch.Start_loc(1);
    batch_end_time = current_batch.Start_loc(end);
    
    % 计算时间差并转换为秒
    batch_duration_samples = batch_end_time - batch_start_time;
    batch_duration_seconds = batch_duration_samples / fs;
    
    % 计算并执行暂停
    pause_duration = batch_duration_seconds * speed_factor;
    pause(pause_duration);
end

hold off;
disp('最终结合版动态绘制完成。');