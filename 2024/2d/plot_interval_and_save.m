%% ==================== 0. 清理与参数设置 ====================
clear; 
clc; 
close all;

% 数据文件路径
filename = 'results\20240822165932_result_yld_3.6e8_5.6e8_window_1024_256_阈值4倍标准差_去零飘_30_80_hann.txt';

% 遍历参数设置
start_base = 3.9e8;       % 起始基准值
end_max = 4e8;          % 最大结束值
step = 0.001e8;            % 步长（每个区间长度）

% 绘图参数
thea = 5000;
point_size = 30;  
num_subplots = 20;  
rows = 4;  % 子图行数
cols = 5;  % 子图列数

% 创建保存图像的文件夹
save_dir = 'step_plots';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

%% ==================== 1. 加载基础数据 ====================
fprintf('--- 正在加载基础数据文件 ---\n');
try
    result_table = readtable(filename);
catch
    error('数据文件 "%s" 加载失败，请检查文件名或路径。', filename);
end

%% ==================== 2. 遍历所有区间并绘图保存 ====================
% 计算总区间数
total_intervals = floor((end_max - start_base) / step);
fprintf('--- 开始遍历 %d 个区间，结果将保存至 %s 文件夹 ---\n', total_intervals, save_dir);

% 循环遍历每个区间
current_start = start_base;
interval_idx = 1;

while current_start < end_max
    current_end = current_start + step;
    
    % 确保最后一个区间不超过最大结束值
    if current_end > end_max
        current_end = end_max;
    end
    
    fprintf('处理区间 %d/%d: [%.3e, %.3e]\n', interval_idx, total_intervals, current_start, current_end);
    
    try
        %% 数据筛选
        logicalIndex = abs(result_table.t123) < 1 & ...
                       abs(result_table.Rcorr) > 0.6 & ...
                       result_table.Start_loc >= current_start & ...
                       result_table.Start_loc < current_end & ...
                       result_table.Elevation < 80 ;
        
        filtered_table = result_table(logicalIndex, :);
        
        if isempty(filtered_table)
            warning('区间 [%.3e, %.3e] 无满足条件的数据，跳过...', current_start, current_end);
            current_start = current_end;
            interval_idx = interval_idx + 1;
            continue;
        end
        
        %% 梯级识别
        filtered_table = sortrows(filtered_table, 'Start_loc');
        time_diffs = diff(filtered_table.Start_loc);
        gap_indices = find(time_diffs > thea);
        
        if length(gap_indices) < 2
            warning('区间 [%.3e, %.3e] 间隙数量不足，跳过...', current_start, current_end);
            current_start = current_end;
            interval_idx = interval_idx + 1;
            continue;
        end
        
        % 确定目标梯级
        gap_duration = diff(gap_indices);
        [~, index] = max(gap_duration);
        step_start_indice = gap_indices(index);
        step_end_indice = gap_indices(index + 1);
        
        current_step_data = filtered_table(step_start_indice+1 : step_end_indice, :);
        num_points = height(current_step_data);
        
        if num_points <= 20
            warning('区间 [%.3e, %.3e] 目标梯级无数据，跳过...', current_start, current_end);
            current_start = current_end;
            interval_idx = interval_idx + 1;
            continue;
        end
        
        %% 数据分份与绘图
        % 定义颜色
        new_stage_colors = lines(num_subplots);
        previous_stage_color = [0.75 0.75 0.75];
        
        % 计算每份点数
        points_per_part = floor(num_points / num_subplots);
        remainder = num_points - points_per_part * num_subplots;
        
        % 计算每份起止索引
        part_indices = cell(num_subplots, 1);
        part_start = 1;
        for p = 1:num_subplots
            part_end = part_start + points_per_part - 1;
            if p <= remainder
                part_end = part_end + 1;
            end
            part_indices{p} = part_start:part_end;
            part_start = part_end + 1;
        end
        
        % 创建图形
        fig = figure('Color', 'w', 'Position', [100 100 1400 900]);
        sgtitle(sprintf('梯级分阶段累积图 (%d - %d)', ...
                        current_step_data.Start_loc(1), current_step_data.Start_loc(end)), ...
                'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
        
        % 绘制子图
        for i = 1:num_subplots
            subplot(rows, cols, i);
            hold on; grid on;
            
            % 绘制前i-1份累积数据（灰色）
            if i > 1  
                prev_indices = [];
                for j = 1:i-1
                    prev_indices = [prev_indices; part_indices{j}'];
                end
                prev_data = current_step_data(prev_indices, :);
                scatter(prev_data.Azimuth, prev_data.Elevation, point_size, previous_stage_color, 'filled');
            end
            
            % 绘制第i份新增数据（彩色）
            current_data = current_step_data(part_indices{i}, :);
            scatter(current_data.Azimuth, current_data.Elevation, point_size, new_stage_colors(i,:), 'filled');
            
            % 坐标轴样式
            if mod(i-1, cols) == 0  % 第一列显示y轴标签
                ylabel('仰角 (°)', 'Color', 'k', 'FontSize', 7);
            end
            if i > rows*(cols-1)  % 最后一行显示x轴标签
                xlabel('方位角 (°)', 'Color', 'k', 'FontSize', 7);
            end
            set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'Box', 'on', ...
                'FontSize', 6, 'TickDir', 'in');
            axis equal;
            hold off;
        end
        
        %% 保存图像
        % 生成文件名（使用整数格式避免科学计数法符号问题）
        start_val = round(current_step_data.Start_loc(1));
        end_val = round(current_step_data.Start_loc(end));
        fig_name = sprintf('%s\\%d_%d.png', save_dir, start_val, end_val);
        print(fig, fig_name, '-dpng', '-r300');  % 高分辨率保存
        close(fig);  % 关闭当前图像释放内存
        
    catch ME
        warning('区间 [%.3e, %.3e] 处理出错: %s，跳过...', current_start, current_end, ME.message);
    end
    
    % 移动到下一个区间
    current_start = current_end;
    interval_idx = interval_idx + 1;
end

fprintf('--- 所有区间处理完成，图像已保存至 %s 文件夹 ---\n', save_dir);