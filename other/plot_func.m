%%  静态图绘制
% --- 1. 数据准备  ---
filename = '20250820151326_1505_result_yld_8e8_11e8_hann_512_128_bandpass_hann_25e6_85e6_.txt';

% 2. 使用 readtable 函数读取数据
%    该函数会自动将第一行作为表头，并根据空格分隔各列
result1 = readtable(filename);
logicalIndex =  abs(result1.t123) < 1  & abs(result1.Rcorr) > 0.4 &  result1.Start_loc < 11e8 & result1.Start_loc > 8e8;
filteredTable1 = result1(logicalIndex, :);


Start_loc = filteredTable1.Start_loc;
colorValues = (Start_loc - min(Start_loc)) / (max(Start_loc) - min(Start_loc)); % 归一化到 [0, 1]

% --- 2. 绘图  ---
% 直接在这里设置figure的深色背景
figure('Color', [0.1 0.1 0.2]); % figure背景设置为深色

% 使用 scatter 绘图，并应用尺寸和透明度优化
scatter(filteredTable1.Azimuth, filteredTable1.Elevation, ...
        2, ... % 尺寸
        colorValues, ...
        'filled', ...
        'MarkerFaceAlpha', 0.6); % 建议加上透明度，深色背景下透明度效果更好

% --- 3. 标签和标题优化 ---
% 设置标题和轴标签的颜色为白色
title('闪电VHF辐射源二维定位图', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'w');
xlabel('方位角 (Azimuth / °)', 'FontSize', 12, 'Color', 'w');
ylabel('仰角 (Elevation / °)', 'FontSize', 12, 'Color', 'w');

% --- 4. 坐标轴和范围设置 ---
% xlim([120, 360]);
% xticks(120:40:360);
% ylim([15, 85]);
% yticks(15:10:85);

% 设置坐标轴的颜色和刻度字体颜色为白色
set(gca, ...
    'FontSize', 11, ...
    'LineWidth', 1.2, ...
    'Color', [0.1 0.1 0.2], ... % Axes背景色和figure背景色保持一致
    'XColor', 'w', ...          % X轴颜色
    'YColor', 'w');             % Y轴颜色

% --- 5. 颜色映射和颜色条优化 ---
colormap('parula'); % 更换为更专业的 colormap
h = colorbar;

% 颜色条标签和刻度颜色为白色
ylabel(h, '归一化发展时间', 'FontSize', 11, 'Color', 'w');
set(h, 'Color', 'w'); % 颜色条的刻度字体颜色

caxis([0, 1]); % 修正颜色范围

% --- 6. 网格和整体风格 ---
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.4, 'Box', 'on'); % 稍微调高GridAlpha，深色背景下看得更清
















signal_length = 3e5;
r_loction = 7e8;
ch1 = read_signal_tdms('20250820151326_1505CH1.tdms',signal_length,r_loction);
figure
plot(ch1)
plot_signal_spectrum(ch1)