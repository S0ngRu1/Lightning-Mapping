%%  静态图绘制
% --- 1. 数据准备  ---
filename = '20250823172542_1505_result_yld_7.4e8_7.9e8_hann_2048_256_bandpass_hann_30e6_90e6.txt';

% 2. 使用 readtable 函数读取数据
%    该函数会自动将第一行作为表头，并根据空格分隔各列
result1 = readtable(filename);
logicalIndex =  abs(result1.t123) < 10  & abs(result1.Rcorr) > 0.1 &  result1.Start_loc < 9e8 & result1.Start_loc > 5.9e8;
filteredTable1 = result1(logicalIndex, :);


Start_loc = filteredTable1.Start_loc;
colorValues = (Start_loc - min(Start_loc)) / (max(Start_loc) - min(Start_loc)); % 归一化到 [0, 1]

% --- 2. 绘图  ---
% 设置figure的浅色背景
figure('Color', [1 1 1]); % figure背景设置为白色

% 使用 scatter 绘图，并应用尺寸和透明度优化
scatter(filteredTable1.Azimuth, filteredTable1.Elevation, ...
        2, ... % 尺寸
        colorValues, ...
        'filled', ...
        'MarkerFaceAlpha', 0.8); % 浅色背景下可适当提高透明度

% --- 3. 标签和标题优化 ---
% 设置标题和轴标签的颜色为深色
title('闪电VHF辐射源二维定位图', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
xlabel('方位角 (Azimuth / °)', 'FontSize', 12, 'Color', 'k');
ylabel('仰角 (Elevation / °)', 'FontSize', 12, 'Color', 'k');

% --- 4. 坐标轴和范围设置 ---
xlim([200, 300]);
xticks(200:20:300);
ylim([5, 85]);
yticks(5:10:85);


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







filename = '标准结果\20250823172542_5453_300000000_200000000_2048_512_8_spectrum-20251021.txt';

% 2. 使用 readtable 函数读取数据
%    该函数会自动将第一行作为表头，并根据空格分隔各列
result1 = readtable(filename);
logicalIndex =  abs(result1.t123) < 0.001  & abs(result1.Rcorr) > 0.001 &  result1.Start_loc < 7.9e8 & result1.Start_loc >7.4e8;
filteredTable1 = result1(logicalIndex, :);


Start_loc = filteredTable1.Start_loc;
colorValues = (Start_loc - min(Start_loc)) / (max(Start_loc) - min(Start_loc)); % 归一化到 [0, 1]

% --- 2. 绘图  ---
% 设置figure的浅色背景
figure('Color', [1 1 1]); % figure背景设置为白色

% 使用 scatter 绘图，并应用尺寸和透明度优化
scatter(filteredTable1.Azimuth, filteredTable1.Elevation, ...
        2, ... % 尺寸
        colorValues, ...
        'filled', ...
        'MarkerFaceAlpha', 0.8); % 浅色背景下可适当提高透明度

% --- 3. 标签和标题优化 ---
% 设置标题和轴标签的颜色为深色
title('闪电VHF辐射源二维定位图', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
xlabel('方位角 (Azimuth / °)', 'FontSize', 12, 'Color', 'k');
ylabel('仰角 (Elevation / °)', 'FontSize', 12, 'Color', 'k');

% --- 4. 坐标轴和范围设置 ---
xlim([200, 300]);
xticks(200:20:300);
ylim([5, 85]);
yticks(5:10:85);

% 设置坐标轴的颜色和刻度字体颜色为深色
set(gca, ...
    'FontSize', 11, ...
    'LineWidth', 1.2, ...
    'Color', [1 1 1], ... % Axes背景色和figure背景色保持一致（白色）
    'XColor', [0.2 0.2 0.2], ... % X轴颜色（深灰色）
    'YColor', [0.2 0.2 0.2]);    % Y轴颜色（深灰色）

% --- 5. 颜色映射和颜色条优化 ---
colormap('parula'); % 保持专业的颜色映射
h = colorbar;

% 颜色条标签和刻度颜色为深色
ylabel(h, '归一化发展时间', 'FontSize', 11, 'Color', 'k');
set(h, 'Color', [0.2 0.2 0.2]); % 颜色条的刻度字体颜色（深灰色）

caxis([0, 1]); % 保持颜色范围

% --- 6. 网格和整体风格 ---
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3, 'Box', 'on'); % 浅色背景下降低网格透明度










signal_length = 3e5;
r_loction = 7e8;
ch1 = read_signal_tdms('20250820151326_1505CH1.tdms',signal_length,r_loction);
figure
plot(ch1)
plot_signal_spectrum(ch1)


% 绘制交替存储的双通道信号
ch1 = read_signal_chj('20250823\software_CardOne_20250823172542.dat', 6e8, 7e8, 2);
x1 = downsample(ch1,50);
plot(x1)

clear
ch4 = read_signal_chj('20250823\software_CardTwo_20250823172542.dat', 2.5e8, 1.3e9+2, 2);
x4 = downsample(ch4,50);
figure
plot(ch4)


clear
yld_ch4 = read_signal_tdms('20250823\20250823172542_5453CH7.tdms',1e8,5.95e8);
plot(yld_ch4)

% 查看两个通道的初始时差
chj_ch1 = read_signal_chj('20250823\software_CardOne_20250823172542.dat', 3e8, 9.6e8, 2);
yld_ch1 = read_signal_tdms('20250823\20250823172542_5453CH3.tdms',3e8,594778000);
figure
plot(chj_ch1)
figure
plot(yld_ch1)



ch1 = read_signal_tdms('20250823\20250823172542_5453CH3.tdms',3e8,5.95e8);
figure
plot(ch1)
ylim([-1000 1000])
plot_signal_spectrum(ch1)



x1 = downsample(processed_ch1_yld,50);
figure
plot(x1)
x2 = downsample(processed_ch2_yld,50);
figure
plot(x2)
x3 = downsample(processed_ch3_yld,50);
figure
plot(x3)



figure
plot(windowed_ch1)
figure
plot(windowed_ch2)
figure
plot(windowed_ch3)



ch1 = read_signal_tdms('20250823\20250823172542_5453CH3.tdms',1024,5.95e8+15063800);
plot(processed_ch1_yld)

plot(processed_ch1_yld(47437900-5120:47437900+5120))
plot_signal_spectrum(processed_ch1_yld(47437900-5120:47437900+5120))