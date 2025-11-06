%% ================== 第 1 部分：加载和处理数据 ==================
% --- 快电场信号 ---
signal_length = 0.8e7;
r_loction_yld = 3.8e8;
ch1_yld = read_signal('..\20240822165932.6610CH4.dat', signal_length, r_loction_yld);
filename = 'results\20240822165932_result_yld_3.6e8_5.6e8_window_1024_256_阈值4倍标准差_去零飘_30_80_hann.txt';
result1 = readtable(filename);
% 信号归一化
y_abs_max = max(abs(ch1_yld));
if y_abs_max ~= 0
    ch1_normalized = ch1_yld / y_abs_max;  % 归一化到[-1,1]
else
    ch1_normalized = ch1_yld;
end
point_size = 15;
sampling_interval_ns = 5;
ns_to_ms = 1e-6; 
time_conversion_factor = sampling_interval_ns * ns_to_ms; % ms/sample

% (修改) X轴时间从 0 开始
x_indices = 0 : signal_length - 1;
time_ms = x_indices * time_conversion_factor;

% --- 1b. 准备VHF辐射源数据 (用于 右Y轴) ---
% result1 已在代码顶部加载或模拟

% --- 1c. 过滤数据并准备绘图 ---
filter_loc_min = r_loction_yld;
filter_loc_max = r_loction_yld+signal_length;
logicalIndex =  abs(result1.t123) < 1  & abs(result1.Rcorr) > 0.6 & ...
                result1.Start_loc < filter_loc_max & result1.Start_loc > filter_loc_min;
filteredTable1 = result1(logicalIndex, :);

% (修改) 将VHF样本位置转换为相对时间 (ms)，使其从 0 开始
vhf_time_ms = (filteredTable1.Start_loc - r_loction_yld) * time_conversion_factor;
vhf_elevation = filteredTable1.Elevation;
% vhf_azimuth 不再需要


%% ================== 第 2 部分：使用 YYAXIS 创建融合图 ==================

figure('Color', [1 1 1], 'Position', [100, 100, 900, 600]);
ax = gca; % 获取当前坐标轴

% --- 2a. 绘制 左Y轴 (快电场) ---
yyaxis(ax, 'left');
plot(ax, time_ms, ch1_normalized, 'Color', 'black', 'LineWidth', 1.5); % 蓝色

% 设置左Y轴
ylabel(ax, '归一化电场强度', 'FontSize', 12, 'FontWeight', 'bold');
set(ax, 'YColor', [0 0.4470 0.7410]); % Y轴标签颜色与线条一致
ylim(ax, [-1.1, 1.1]); % 归一化范围


% --- 2b. 绘制 右Y轴 (VHF仰角) ---
yyaxis(ax, 'right');
% 颜色参数 (第5个) 改为 vhf_time_ms
scatter(ax, vhf_time_ms, vhf_elevation, ...
        point_size, ...           % 尺寸
        vhf_time_ms, ...  % 颜色数据 (时间)
        'filled', ...
        'MarkerFaceAlpha', 0.7);

% 设置右Y轴
ylabel(ax, '仰角 (Elevation / °)', 'FontSize', 12, 'FontWeight', 'bold');
set(ax, 'YColor', [0.2 0.2 0.2]); % 右Y轴使用深灰色
ylim(ax, [0, 85]); % 您原始的仰角范围
yticks(ax, 0:10:85);


% --- 2c. 格式化 共享X轴 (时间) ---
xlabel(ax, '时间 (ms)', 'FontSize', 12, 'FontWeight', 'bold');

% (修改) X轴范围从 0 开始
time_min_ms = 0;
time_max_ms = (signal_length - 1) * time_conversion_factor;
xlim(ax, [time_min_ms, time_max_ms]);

% --- 2d. 格式化 颜色条 (时间) ---
h = colorbar(ax);
% 颜色条标签改为 "时间 (ms)"
ylabel(h, '时间 (ms)', 'FontSize', 11);
colormap(ax, 'parula'); % 或 'jet', 'hsv'
% 颜色范围对应于X轴的时间范围
clim(ax, [time_min_ms, time_max_ms]); 

% --- 2e. 最终修饰 ---
title('闪电快电场与VHF辐射源时序图', 'FontSize', 16, 'FontWeight', 'bold');
grid(ax, 'on');
set(ax, 'FontSize', 11, 'LineWidth', 1.2, 'Box', 'on', 'GridAlpha', 0.3, 'GridLineStyle', '--');

% (更新状态描述
disp('yyaxis 融合图绘制完成。');
disp('X轴 (下): 时间 (ms)');
disp('Y轴 (左): 归一化电场 (蓝色)');
disp('Y轴 (右): 仰角 (黑色)');
disp('散点颜色: 时间 (ms) (Parula色图)');