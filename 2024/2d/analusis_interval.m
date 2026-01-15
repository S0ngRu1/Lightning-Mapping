%% --- 1. 初始化和数据读取 ---
clear; clc; close all;

% 将您的数据保存到名为 'data.txt' 的文件中
filename = 'results\result_yld_3.8e8_4e8_window_512_128_去零飘_滤波_加窗_阈值15_30_80.txt';

% 使用 readtable 读取数据，它会自动处理表头和空格分隔符
try
    dataTable = readtable(filename);
catch ME
    error('无法读取文件 "%s"。请确保文件存在于当前目录，并且格式正确。', filename);
end

logicalIndex =  abs(dataTable.t123) < 1  & abs(dataTable.Rcorr) > 0.3 &  dataTable.Start_loc < 4e8 & dataTable.Start_loc > 3.8e8;
filteredTable1 = dataTable(logicalIndex, :);


%% --- 2. 定义时间区间并筛选数据 ---
interval1_start = 3.8e8;
interval1_end = 3.88e8;
interval2_start = 3.88e8;
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

%% --- 3. 分析正负先导的发展特征 ---
for leader_case = 1:2

    % --- 2.1 根据循环次数，选择要处理的数据和参数 ---
    if leader_case == 1
        % --- 第一次循环：处理正先导 ---
        current_table = table_interval1;
        leader_name = '正先导';
        time_diffs = diff(current_table.Start_loc);
        figure;
        h = histogram(time_diffs);
        set(gca, 'XScale', 'log');
        title('正先导的时间差分布 (用于确定阈值)');
        xlabel('Time Difference (sampling points) - Log Scale');
        ylabel('Frequency (Count)');
        grid on;

        % 正先导是连续发展的，我们将其视为一个完整的阶段
        if ~isempty(current_table)
            step_start_indices = 1;
            step_end_indices = height(current_table);
        else
            step_start_indices = [];
            step_end_indices = [];
        end

    else
        % --- 第二次循环：处理负先导 ---
        current_table = table_interval2;
        leader_name = '负先导';

        % 对负先导进行详细的阶段划分分析
        if height(current_table) > 1
            time_diffs = diff(current_table.Start_loc);

            % 【注意】直方图仅为负先导生成
            figure;
            h = histogram(time_diffs);
            set(gca, 'XScale', 'log');
            title('负先导的时间差分布 (用于确定阈值)');
            xlabel('Time Difference (sampling points) - Log Scale');
            ylabel('Frequency (Count)');
            grid on;

            gap_threshold = 39000;
            gap_indices = find(time_diffs > gap_threshold);
            step_start_indices = [1; gap_indices + 1];
            step_end_indices = [gap_indices; height(current_table)];
        elseif ~isempty(current_table)
            step_start_indices = 1;
            step_end_indices = 1;
        else
            step_start_indices = [];
            step_end_indices = [];
        end
    end

%     % 如果当前数据表为空，则跳过本次循环
%     if isempty(current_table)
%         fprintf('\n%s 的数据为空，跳过绘制。\n', leader_name);
%         continue;
%     end
% 
%     %% --- 2.2 通用的阶段分配和绘图逻辑 ---
%     % (这部分代码现在是通用的，无需再复制)
% 
%     % a. 为每个定位点分配一个阶段ID
%     num_steps = numel(step_start_indices);
%     current_table.StepID = zeros(height(current_table), 1);
%     for i = 1:num_steps
%         start_idx = step_start_indices(i);
%         end_idx = step_end_indices(i);
%         current_table.StepID(start_idx:end_idx) = i;
%     end
% 
%     % b. 为每个阶段ID创建一种独特的颜色
%     step_colors = lines(num_steps);
% 
%     % c. 参数定义
%     fs = 200e6;
%     speed_factor = 500;
%     plot_batch_size = 20;
% 
%     % d. 图形初始化 (每次循环都会创建一个新的figure窗口)
%     figure;
%     hold on;
%     grid on;
%     box on;
%     xlabel('方位角 (Azimuth)');
%     xlim([120, 200]);
%     xticks(120:20:200);
%     ylabel('仰角 (Elevation)');
%     ylim([10, 70]);
%     yticks(10:10:70);
%     title(sprintf('%s 发展阶段动态呈现', leader_name)); % 使用变量设置标题
% 
%     % e. 创建图例
%     legend_handles = gobjects(num_steps, 1);
%     legend_entries = cell(num_steps, 1);
%     for i = 1:num_steps
%         legend_handles(i) = scatter(NaN, NaN, 50, step_colors(i,:), 'filled');
%         legend_entries{i} = sprintf('%s 阶段 %d', leader_name, i);
%     end
%     legend(legend_handles, legend_entries, 'Location', 'northwest');
% 
%     % f. 动态绘制循环
%     for i = 1:plot_batch_size:height(current_table)
%         start_idx = i;
%         end_idx = min(i + plot_batch_size - 1, height(current_table));
%         current_batch = current_table(start_idx:end_idx, :);
% 
%         if isempty(current_batch)
%             continue;
%         end
% 
%         batch_step_ids = current_batch.StepID;
%         valid_ids = batch_step_ids > 0;
%         if any(valid_ids)
%             batch_colors = step_colors(batch_step_ids(valid_ids), :);
%             scatter(current_batch.Azimuth(valid_ids), current_batch.Elevation(valid_ids), 5, batch_colors, 'filled', 'HandleVisibility', 'off');
%         end
% 
%         drawnow;
% 
%         batch_start_time = current_batch.Start_loc(1);
%         batch_end_time = current_batch.Start_loc(end);
%         batch_duration_samples = batch_end_time - batch_start_time;
%         batch_duration_seconds = batch_duration_samples / fs;
%         pause_duration = batch_duration_seconds * speed_factor;
%         pause(pause_duration);
%     end
% 
%     hold off;
%     fprintf('\n%s 动态绘制完成。\n', leader_name);

end % 循环结束