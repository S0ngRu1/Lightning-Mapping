%% 自动计算 20 个有效区间 (Intervals Generation - Target 20)
clc; clear; close all;

% -------------------------- 1. 数据读取与预处理 --------------------------
filename = 'results\20240822165932_result_yld_3.6e8_5.6e8_window_1024_256_阈值4倍标准差_去零飘_30_80_hann.txt';
Start_loc_Base = 3.709e8; % 基准值
if ~isfile(filename), error('文件不存在'); end
base_value = Start_loc_Base;
point_size = 15;
stepSize = 1;
% 读取数据
result1 = readtable(filename);

% 筛选数据
logicalIndex =  abs(result1.t123) < 0.5  & ...
                abs(result1.Rcorr) > 0.5 & ...
                result1.Start_loc < Start_loc_Base + 5e5 & ...
                result1.Start_loc > Start_loc_Base & ...
                result1.Elevation < 20 & ...
                result1.Azimuth > 180;
filteredTable1 = result1(logicalIndex, :);
% 提取 offset 数据
valid_start_locs = filteredTable1.Start_loc;
offset = valid_start_locs - Start_loc_Base;

if isempty(offset)
    error('筛选后数据为空，无法计算区间。');
end

% 计算相邻点之间的差值 (Gaps)
diff_start_loc = diff(offset);

% -------------------------- 2. 动态计算区间 (核心算法) --------------------------
% 目标：找到 20 个区间，且每个区间长度 (End - Start) > 2

target_count = 20;       % 目标区间数量
min_interval_len = 2;    % 最小区间长度限制

% 对所有间隙从大到小排序，获取排序后的索引
[~, sortIdx] = sort(diff_start_loc, 'descend');

found_intervals = [];
final_intervals = [];

max_possible_cuts = length(diff_start_loc);
start_cuts = target_count - 1; 

for k = start_cuts : max_possible_cuts
    % 1. 取出当前最大的 k 个间隙的索引
    % 必须重新按索引大小排序(sort)，以保证时间轴的先后顺序
    current_gap_indices = sort(sortIdx(1:k));
    
    % 2. 根据这些断点构建临时区间
    % 中间的 End 和 Start
    middle_ends = offset(current_gap_indices) + 1; 
    middle_starts = offset(current_gap_indices + 1) - 1;
    
    % 首尾
    first_start = offset(1) - 1;
    last_end = offset(end) + 1;
    
    % 拼接
    temp_starts = [first_start; middle_starts];
    temp_ends   = [middle_ends; last_end];
    temp_intervals = [temp_starts, temp_ends];
    
    % 3. 筛选：检查长度是否满足条件 > 2
    % 计算长度：End - Start
    lens = temp_intervals(:, 2) - temp_intervals(:, 1);
    valid_mask = lens > min_interval_len;
    
    valid_intervals = temp_intervals(valid_mask, :);
    
    % 4. 检查数量是否达标
    if size(valid_intervals, 1) >= target_count
        % 如果有效区间数量达到或超过目标，取前 target_count 个（按时间顺序）
        final_intervals = valid_intervals(1:target_count, :);
        disp(['成功找到 ', num2str(target_count), ' 个有效区间！(使用了 ', num2str(k), ' 个切分点)']);
        break;
    end
end

if isempty(final_intervals) || size(final_intervals, 1) < target_count
    warning(['未能找到足够的区间。仅找到了 ', num2str(size(valid_intervals, 1)), ' 个符合条件的区间。']);
    final_intervals = valid_intervals; % 返回尽力找到的结果
end
intervals = final_intervals;
disp(intervals);

%% 精细化绘制：带全局统计标题的闪电发展过程

% intervals = [
%         4996        7548
%        23780       25358
%        44612       47444
%        64419       65623
%        67130       70267
%        73757       78134
%        81112       82340
%       100844      104152
%       111151      112261
%       120178      122024
%       126467      128155
%       138195      139745
%       147062      149721
%       154634      156109
%       162719      169446
%       179042      182158
%       184527      188775
% ];
n = size(intervals, 1); 

% -------------------------- 2. 统计计算 --------------------------
% A. 计算所有区间的持续时间 (Duration)
raw_durs = intervals(:, 2) - intervals(:, 1);
us_durs = raw_durs * 5 / 1e3; 
% B. 计算所有区间的梯级间隔 (Gap)
if n > 1
    raw_gaps = intervals(2:end, 1) - intervals(1:end-1, 2);
    us_gaps = raw_gaps * 5 / 1e3; 
else
    us_gaps = 0; 
end
% C. 生成标题
titleStr_Dur = sprintf('持续时间：%.1fμs (%.1f-%.1fμs)', mean(us_durs), min(us_durs), max(us_durs));
if n > 1
    titleStr_Gap = sprintf('梯级间隔：%.1fμs (%.1f-%.1fμs)', mean(us_gaps), min(us_gaps), max(us_gaps));
else
    titleStr_Gap = '梯级间隔：N/A';
end
fullTitle = [titleStr_Dur, ',  ', titleStr_Gap];

% -------------------------- 3. 视频录制初始化 --------------------------
videoName = 'Lightning_Statistics_FixedBackground.mp4'; 
v = VideoWriter(videoName, 'MPEG-4');   
v.FrameRate = 2; 
v.Quality = 95;  
open(v);         

% -------------------------- 4. 数据读取 --------------------------
if ~isfile(filename), error('文件不存在'); end
result1 = readtable(filename);
baseData = filteredTable1;
az_min = min(filteredTable1.Azimuth);
az_max = max(filteredTable1.Azimuth);
el_min = min(filteredTable1.Elevation);
el_max = max(filteredTable1.Elevation);

% -------------------------- 5. 绘图窗口初始化 --------------------------
hFig = figure('Name', '闪电区间统计', 'Color', 'w', 'Visible', 'on');
set(hFig, 'Position', [50, 50, 1600, 900]); 
colors = lines(n); 
sgtitle(fullTitle, 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');

% ========================== 循环处理每个子图 ==========================
for i = 1:n
    ax = subplot(4, 5, i);
    hold(ax, 'on');
    
    % --- 坐标轴设置 ---

    xlim([az_min, az_max]);  % 固定x轴范围（动态绘图不变化）
    xticks('auto');          % 让Matlab在固定范围内自动生成规整刻度
    ylim([el_min, el_max]);  % 固定y轴范围
    yticks('auto');  
    grid(ax, 'on');
    set(ax, 'FontSize', 8, 'LineWidth', 1, 'Box', 'on', ...
        'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2], ...
        'GridLineStyle', '--', 'GridAlpha', 0.3);
    
    % --- 计算当前子图时间范围 ---
    curr_start_loc = base_value + intervals(i, 1);
    curr_end_loc   = base_value + intervals(i, 2);
    
    % --- 标签 ---
    this_dur = (intervals(i, 2) - intervals(i, 1)) * 5 / 1e3;
    if i == 1
        xlabelStr = ['Dur: ', num2str(this_dur, '%.1f'), ' μs'];
    else
        this_gap = (intervals(i, 1) - intervals(i-1, 2)) * 5 / 1e3;
        xlabelStr = {['Dur: ', num2str(this_dur, '%.1f'), ' μs']; ...
                     ['Gap: ', num2str(this_gap, '%.1f'), ' μs']};
    end
    
    title(ax, ['区间 ', num2str(i)], 'FontWeight', 'bold', 'FontSize', 10);
    xlabel(ax, xlabelStr, 'FontSize', 9, 'Color', 'b'); 
    ylabel(ax, 'El (°)');

    if i > 1
        % 初始化一个全假的掩码
        historyMask = false(height(baseData), 1);
        
        % 循环累加前 i-1 个区间的有效索引
        for k = 1:i-1
            k_start = base_value + intervals(k, 1);
            k_end   = base_value + intervals(k, 2);
            % 逻辑或运算：只要点属于之前任何一个区间，就标记为 True
            historyMask = historyMask | (baseData.Start_loc >= k_start & baseData.Start_loc <= k_end);
        end
        
        % 提取历史数据
        prevData = baseData(historyMask, :);
        
        if ~isempty(prevData)
            scatter(ax, prevData.Azimuth, prevData.Elevation, point_size, [0.6 0.6 0.6], ...
                'filled', 'MarkerFaceAlpha', 0.3);
        end
    end
    
    drawnow;
    writeVideo(v, getframe(hFig));
    
    % --- 绘制当前区间动态数据 ---
    newIndex = baseData.Start_loc >= curr_start_loc & baseData.Start_loc <= curr_end_loc;
    newData = baseData(newIndex, :);
    numPoints = height(newData);
    
    if numPoints > 0
        newData = sortrows(newData, 'Start_loc');
        allAz = newData.Azimuth;
        allEl = newData.Elevation;
        
        hLine = plot(ax, NaN, NaN, '-', 'Color', [colors(i,:) 0.4], 'LineWidth', 0.5);
        hScatter = scatter(ax, NaN, NaN, point_size, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.8);
        hHead = plot(ax, NaN, NaN, 'r+', 'MarkerSize', point_size, 'LineWidth', 1.5);
        
        accX = [];
        accY = [];
        
        for j = 1:stepSize:numPoints
            idxStart = j;
            idxEnd = min(j + stepSize - 1, numPoints);
            
            chunkX = allAz(idxStart:idxEnd);
            chunkY = allEl(idxStart:idxEnd);
            
            accX = [accX; chunkX];
            accY = [accY; chunkY];
            
            set(hLine, 'XData', accX, 'YData', accY);
            set(hScatter, 'XData', accX, 'YData', accY);
            set(hHead, 'XData', chunkX(end), 'YData', chunkY(end));
            
            % Start/End 标记
            if j == 1
                scatter(ax, chunkX(1), chunkY(1), point_size*2, 'g', '^', 'filled', 'MarkerEdgeColor', 'k');
            end
            if idxEnd == numPoints
                scatter(ax, chunkX(end), chunkY(end), point_size*2, 'r', 's', 'filled', 'MarkerEdgeColor', 'k');
            end
            
            drawnow;
            writeVideo(v, getframe(hFig));
        end
        set(hHead, 'Visible', 'off');
    else
        writeVideo(v, getframe(hFig));
    end
    
    hold(ax, 'off');
end

% -------------------------- 结束处理 --------------------------
close(v);
disp('视频录制完成！');