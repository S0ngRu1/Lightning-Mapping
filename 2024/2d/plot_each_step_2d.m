%% 根据每一个梯级的起始位置绘制每一个梯级的二维结果
% -------------------------- 参数设置 --------------------------
filename = 'results\20240822165932_result_yld_3.65e8_5e8_window_256_64_阈值4倍标准差_去零飘_30_80_hann.txt';
base_value = 4.176e8;  % 基准值（采样位置起点）
start_positions = [0,884, 1503, 31402, 33394, 35094, 37558];  % 起始位置列表
n = length(start_positions);  % 总区间数，最终绘制n-1个子图

% 预计算每个区间的独立持续时间（diff(start_positions)→相邻起点差，*5/1e3转μs）
durations = diff(start_positions) * 5 / 1e3;  

% -------------------------- 数据读取与预处理 --------------------------
result1 = readtable(filename);
% 先筛选出所有满足基础条件的数据（t123和Rcorr的条件）
baseIndex = abs(result1.t123) < 1 & abs(result1.Rcorr) > 0.6;
baseData = result1(baseIndex, :);

% 预计算每个区间的上下限（采样位置）
bounds = base_value + start_positions;  % 所有区间的分界点（含起点和终点）

% -------------------------- 绘制累积叠加子图 --------------------------
figure; % 放大窗口避免重叠
% 定义每个新增区间的颜色（区分累积的不同阶段）
colors = lines(n-1);  % 生成n-1种区分度较高的颜色

for i = 1:n-1
    % 当前子图包含第1到第i个区间的数据（累积叠加）
    currentIndex = baseData.Start_loc > bounds(1) & baseData.Start_loc < bounds(i+1);
    currentData = baseData(currentIndex, :);
    
    % 子图布局（2行5列，适配10个子图）
    ax = subplot(2,3, i); % 获取子图句柄，方便后续设置
    hold on; % 统一hold on，避免重复切换
    
    % 分阶段绘制：先画之前累积的区间（灰色，半透明）
    if i > 1
        prevIndex = baseData.Start_loc > bounds(1) & baseData.Start_loc < bounds(i);
        prevData = baseData(prevIndex, :);
        scatter(prevData.Azimuth, prevData.Elevation, 10, [0.5 0.5 0.5], ...
            'filled', 'MarkerFaceAlpha', 0.3);
    end
    
    % 绘制第i个新增区间的数据（用对应颜色）
    newIndex = baseData.Start_loc > bounds(i) & baseData.Start_loc < bounds(i+1);
    newData = baseData(newIndex, :);
    scatter(newData.Azimuth, newData.Elevation, 10, colors(i,:), ...
        'filled', 'MarkerFaceAlpha', 0.8);
    
    % -------------------------- 子图标题与坐标轴设置 --------------------------
    % 子图标题：当前区间范围（[bounds(i), bounds(i+1))）
    title(['区间: [', num2str(bounds(i)), ', ', num2str(bounds(i+1)), ')'], ...
        'FontSize', 10, 'FontWeight', 'bold');
    
    % 横坐标标题：当前区间的独立持续时间（diff(start_positions)*5/1e3）
    xlabel(['持续时间: ', num2str(durations(i), '%.2f'), ' μs'], ...
        'FontSize', 9, 'Color', 'k');
    
    % 纵坐标保留仰角标签
    ylabel('仰角 (°)', 'FontSize', 9, 'Color', 'k');
    
    % 坐标轴范围与刻度
    xlim([192, 204]);
    xticks(192:2:204); % 保留方位角刻度数值（数据本身是方位角）
    ylim([54, 70]);
    yticks(54:2:70);
    
    % 坐标轴样式优化
    set(ax, ...
        'FontSize', 8, ...
        'LineWidth', 1, ...
        'Color', [1 1 1], ...
        'XColor', [0.2 0.2 0.2], ...
        'YColor', [0.2 0.2 0.2], ...
        'Position', get(ax, 'Position') + [0, 0.02, 0, -0.02]); % 微调位置
    
    % 网格设置
    grid on;
    set(ax, 'GridLineStyle', '--', 'GridAlpha', 0.3, 'Box', 'on');
    
    hold off;
end
