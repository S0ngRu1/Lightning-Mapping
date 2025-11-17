%% 根据每一个梯级的起始位置绘制每一个梯级的二维结果
% -------------------------- 参数设置 --------------------------
filename = 'results\20240822165932_result_yld_3.65e8_4.05e8_window_256_64_阈值4倍标准差_去零飘_30_80_hann.txt';
base_value = 365756907;  % 基准值
start_positions = [955, 1951, 2350, 3855, 4540, 5623, 8033, 9518, 11995, 13433];  % 起始位置列表
n = length(start_positions);  % 总区间数，最终绘制n-1个子图

% -------------------------- 数据读取与预处理 --------------------------
result1 = readtable(filename);
% 先筛选出所有满足基础条件的数据（t123和Rcorr的条件）
baseIndex = abs(result1.t123) < 1 & abs(result1.Rcorr) > 0.6;
baseData = result1(baseIndex, :);

% 预计算每个区间的上下限（用于后续筛选）
bounds = base_value + start_positions;  % 所有区间的分界点（含起点和终点）

% -------------------------- 绘制累积叠加子图 --------------------------
figure('Color', [1 1 1], 'Position', [100 100 1200 800]);
% 定义每个新增区间的颜色（区分累积的不同阶段）
colors = lines(n-1);  % 生成n-1种区分度较高的颜色

for i = 1:n-1
    % 当前子图包含第1到第i个区间的数据（累积叠加）
    % 区间范围：[bounds(1), bounds(i+1))
    currentIndex = baseData.Start_loc > bounds(1) & baseData.Start_loc < bounds(i+1);
    currentData = baseData(currentIndex, :);
    
    % 子图布局（2行3列，可根据数量调整）
    subplot(3, 3, i);
    
    % 分阶段绘制：先画之前累积的区间（灰色），再画新增区间（彩色）
    if i > 1
        % 绘制1到i-1区间的累积数据（灰色，半透明）
        prevIndex = baseData.Start_loc > bounds(1) & baseData.Start_loc < bounds(i);
        prevData = baseData(prevIndex, :);
        scatter(prevData.Azimuth, prevData.Elevation, 20, [0.5 0.5 0.5], ...
            'filled', 'MarkerFaceAlpha', 0.3);
        hold on;  % 保持当前图，用于叠加新数据
    end
    
    % 绘制第i个新增区间的数据（用对应颜色）
    newIndex = baseData.Start_loc > bounds(i) & baseData.Start_loc < bounds(i+1);
    newData = baseData(newIndex, :);
    scatter(newData.Azimuth, newData.Elevation, 20, colors(i,:), ...
        'filled', 'MarkerFaceAlpha', 0.8);
    
    % 子图标题（显示累积区间范围）
    title(['区间: [', num2str(bounds(1)), ', ', num2str(bounds(i+1)), ')'], ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
    
    % 坐标轴设置
    xlabel('方位角 (°)', 'FontSize', 9, 'Color', 'k');
    ylabel('仰角 (°)', 'FontSize', 9, 'Color', 'k');
    xlim([178, 186]);
    xticks(178:1:186);
    ylim([45, 50]);
    yticks(45:0.5:50);
    % 坐标轴样式
    set(gca, ...
        'FontSize', 8, ...
        'LineWidth', 1, ...
        'Color', [1 1 1], ...
        'XColor', [0.2 0.2 0.2], ...
        'YColor', [0.2 0.2 0.2]);
    
    % 网格设置
    grid on;
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3, 'Box', 'on');
    hold off;  % 结束当前子图的叠加
end

% 整体标题与图例
sgtitle('闪电VHF辐射源累积发展过程（按区间叠加）', 'FontSize', 16, 'FontWeight', 'bold');

% 添加图例（说明新增区间的颜色对应关系）
legendAx = axes('Position', [0.05, 0.02, 0.9, 0.05], 'Visible', 'off');
for i = 1:n-1
    scatter(legendAx, i, 0.5, 20, colors(i,:), 'filled', 'DisplayName', ['第', num2str(i), '区间']);
    hold on;
end
legend('Location', 'eastoutside', 'FontSize', 8, 'Orientation', 'horizontal');