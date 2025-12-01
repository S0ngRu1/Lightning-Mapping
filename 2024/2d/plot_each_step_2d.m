%% 根据每一个梯级的起始位置绘制每一个梯级的二维结果
% -------------------------- 参数设置 --------------------------
filename = 'results\20240822165932_result_yld_3.65e8_5e8_window_256_64_阈值4倍标准差_去零飘_30_80_hann.txt';
base_value = 4.176e8;  % 基准值（采样位置起点）
start_positions = [0,884, 1503, 31402, 33394, 35094, 37558];  % 起始位置列表
n = length(start_positions);  % 总区间数，最终绘制n-1个子图

% 预计算每个区间的独立持续时间
durations = diff(start_positions) * 5 / 1e3;  

% --- 动态绘制速度控制 ---
% stepSize 控制每次刷新屏幕时新增点的数量。
% stepSize = 1 会最流畅但最慢；设置大一些（如 20-50）可以加快录制速度，但颗粒感会增加。
stepSize = 1; 

% -------------------------- 视频录制初始化 --------------------------
videoName = 'Dynamic_Scatter_Growth.mp4'; % 视频保存文件名
v = VideoWriter(videoName, 'MPEG-4');   % 创建视频对象
% 帧率决定了视频播放的速度。因为我们是在内循环中抓帧，实际生成的帧数非常多。
% 设置较高的帧率可以让播放看起来更顺滑。
v.FrameRate = 3; 
v.Quality = 95;  % 视频质量 (0-100)
open(v);         % 打开视频文件准备写入

disp(['正在创建动态视频: ', videoName, ' ...']);
disp('录制过程可能需要几分钟，请耐心等待，不要移动图形窗口。');

% -------------------------- 数据读取与预处理 --------------------------
result1 = readtable(filename);
% 筛选数据
baseIndex = abs(result1.t123) < 1 & abs(result1.Rcorr) > 0.6;
baseData = result1(baseIndex, :);
% 预计算区间上下限
bounds = base_value + start_positions;

% -------------------------- 初始化图形窗口 --------------------------
% 创建一个固定大小的窗口，保证视频分辨率一致
hFig = figure('Name', '动态散点生长过程', 'Color', 'w', 'Visible', 'on');
% 设置窗口位置和大小 [left, bottom, width, height]，例如 1600x900 分辨率
set(hFig, 'Position', [50, 50, 1600, 900]); 

% 定义颜色
colors = lines(n-1); 

% ========================== 外层循环：遍历子图 ==========================
for i = 1:n-1
    % --- 1. 激活并设置当前子图坐标轴 ---
    ax = subplot(2, 3, i);
    hold(ax, 'on');
    
    % 预先设置好坐标轴范围和标题，避免动画过程中坐标轴跳动
    xlim(ax, [192, 204]); xticks(ax, 192:2:204);
    ylim(ax, [54, 70]);   yticks(ax, 54:2:70);
    grid(ax, 'on');
    set(ax, 'FontSize', 8, 'LineWidth', 1, 'Box', 'on', ...
        'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2], ...
        'GridLineStyle', '--', 'GridAlpha', 0.3);
    
    title(ax, ['区间: [', num2str(bounds(i)), ', ', num2str(bounds(i+1)), ')'], ...
        'FontSize', 10, 'FontWeight', 'bold');
    xlabel(ax, ['持续时间: ', num2str(durations(i), '%.2f'), ' μs'], ...
        'FontSize', 9);
    ylabel(ax, '仰角 (°)', 'FontSize', 9);

    % --- 2. 绘制背景（之前区间的累积数据，灰色） ---
    if i > 1
        prevIndex = baseData.Start_loc > bounds(1) & baseData.Start_loc < bounds(i);
        prevData = baseData(prevIndex, :);
        if ~isempty(prevData)
            % 一次性绘制所有背景点，不需要动态
            scatter(ax, prevData.Azimuth, prevData.Elevation, 10, [0.6 0.6 0.6], ...
                'filled', 'MarkerFaceAlpha', 0.3);
        end
    end
    
    % 在开始动态绘制当前区间前，先抓取一帧作为初始状态
    drawnow;
    writeVideo(v, getframe(hFig));
    
    % --- 3. 准备当前区间的动态数据 ---
    newIndex = baseData.Start_loc > bounds(i) & baseData.Start_loc < bounds(i+1);
    newData = baseData(newIndex, :);
    numPoints = height(newData);
    
    if numPoints > 0
        % 初始化一个空的散点句柄，用于后续动态更新
        % 使用 NaN 初始化，这样一开始图上什么都没有
        hDynamicScatter = scatter(ax, NaN, NaN, 10, colors(i,:), ...
            'filled', 'MarkerFaceAlpha', 0.8);
        
        % 提取当前需要绘制的所有数据
        allAzimuth = newData.Azimuth;
        allElevation = newData.Elevation;
        
        currentXData = [];
        currentYData = [];
        
        disp(['  -> 正在动态绘制第 ', num2str(i), ' 个子图，共 ', num2str(numPoints), ' 个点...']);
        
        % ========================== 内层循环：动态绘制点 ==========================
        % 使用 stepSize 进行批处理更新，提高效率
        for j = 1:stepSize:numPoints
            % 确定当前批次的结束索引
            endIdx = min(j + stepSize - 1, numPoints);
            
            % 追加新的数据点
            currentXData = [currentXData; allAzimuth(j:endIdx)];
            currentYData = [currentYData; allElevation(j:endIdx)];
            
            % --- 核心：更新散点图的数据源 ---
            set(hDynamicScatter, 'XData', currentXData, 'YData', currentYData);
            
            % --- 核心：刷新屏幕并录制 ---
            drawnow; % 强制刷新显示
            frame = getframe(hFig); % 捕获整个窗口
            writeVideo(v, frame);   % 写入视频
        end
    end
    
    hold(ax, 'off');
    disp(['完成第 ', num2str(i), ' 个子图绘制。']);
end

% -------------------------- 结束处理 --------------------------
close(v); % 关闭视频流，完成文件保存
disp('------------------------------------------------');
disp(['视频录制完成！文件已保存为: ', videoName]);