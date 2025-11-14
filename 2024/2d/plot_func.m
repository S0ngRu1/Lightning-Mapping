%%  静态图绘制
% --- 1. 数据准备  ---
filename = 'results\20240822165932_result_yld_3.65e8_4.05e8_window_256_64_阈值4倍标准差_去零飘_30_80_hann.txt';

% 2. 使用 readtable 函数读取数据
%    该函数会自动将第一行作为表头，并根据空格分隔各列
result1 = readtable(filename);
logicalIndex =  abs(result1.t123) < 1  & abs(result1.Rcorr) > 0.6 &  result1.Start_loc < 365756907+13731 & result1.Start_loc > 365756907;
filteredTable1 = result1(logicalIndex, :);


Start_loc = filteredTable1.Start_loc;
colorValues = (Start_loc - min(Start_loc)) / (max(Start_loc) - min(Start_loc)); % 归一化到 [0, 1]

% --- 2. 绘图  ---
% 设置figure的浅色背景
figure('Color', [1 1 1]); % figure背景设置为白色

% 使用 scatter 绘图，并应用尺寸和透明度优化
scatter(filteredTable1.Azimuth, filteredTable1.Elevation, ...
        20, ... % 尺寸
        colorValues, ...
        'filled', ...
        'MarkerFaceAlpha', 0.8); % 浅色背景下可适当提高透明度

% --- 3. 标签和标题优化 ---
% 设置标题和轴标签的颜色为深色
title('闪电VHF辐射源二维定位图', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
xlabel('方位角 (Azimuth / °)', 'FontSize', 12, 'Color', 'k');
ylabel('仰角 (Elevation / °)', 'FontSize', 12, 'Color', 'k');

% --- 4. 坐标轴和范围设置 ---
xlim([120, 220]);
xticks(120:20:220);
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




%% 按比例的“实时回放”
filename = 'results\20240822165932_result_yld_3.65e8_4.05e8_window_256_64_阈值4倍标准差_去零飘_30_80_hann.txt';

% 2. 使用 readtable 函数读取数据
%    该函数会自动将第一行作为表头，并根据空格分隔各列
result1 = readtable(filename);
logicalIndex =  abs(result1.t123) < 1  & abs(result1.Rcorr) > 0.6 &  result1.Start_loc < 365756907+13731 & result1.Start_loc > 365756907;
filteredTable1 = result1(logicalIndex, :);


Fs = 200e6; % 采样率 
target_viz_duration = 40; % (秒) 播放总时长
min_viz_pause = 5e-9; % (秒) 最小暂停时间

% 2. 计算真实时间和颜色
Start_loc = filteredTable1.Start_loc;
Start_loc_min = min(Start_loc);
Start_loc_max = max(Start_loc);

all_event_times_real = (Start_loc - Start_loc_min) / Fs;
all_colorValues = (Start_loc - Start_loc_min) / (Start_loc_max - Start_loc_min);
T_duration_real = max(all_event_times_real); % 闪电事件实际总时长

disp(['闪电事件实际持续时间: ', num2str(T_duration_real * 1e6), ' 微秒 (us)']);
disp(['将缩放到 ', num2str(target_viz_duration), ' 秒进行可视化回放...']);

% 3. 按时间顺序，对所有数据点进行排序
[sorted_times, sort_idx] = sort(all_event_times_real);
sorted_az = filteredTable1.Azimuth(sort_idx);
sorted_el = filteredTable1.Elevation(sort_idx);
sorted_colors = all_colorValues(sort_idx);

% 4. 计算点与点之间的真实时间间隔
% diff([0; V]) 会计算 t1-0, t2-t1, t3-t2 ...
time_gaps_real = diff([0; sorted_times]);

% 5. 计算缩放后的可视化暂停时间
% 缩放因子 = 目标总时长 / 真实总时长
scale_factor = target_viz_duration / T_duration_real;
viz_pauses = time_gaps_real * scale_factor;

% 6. 创建图形窗口
figure('Position', [100, 100, 800, 600]);
hold on;
grid on;
xlabel('方位角');
xlim([178, 186]);
xticks(178:1:186);
ylabel('仰角');
ylim([45, 50]);
yticks(45:0.5:50);
title('目标点动态呈现 (按比例实时回放)');
colormap('hsv');
h_bar = colorbar;
ylabel(h_bar, '归一化起始位置 (代表时间)');
caxis([0, 1]);

% 7. 逐点绘制，并按比例暂停
N = length(sorted_times);
disp(['开始回放... 总共 ', num2str(N), ' 个点。']);
tic; % 开始计时

for i = 1:N
    % (1) 获取当前点的暂停时间
    current_pause = viz_pauses(i);
    
    % (2) 执行暂停
    % 只有当暂停时间大于最小阈值时才真正暂停
    if current_pause > min_viz_pause
        pause(current_pause);
    end
    
    % (3) 绘制当前点
    scatter(sorted_az(i), sorted_el(i), 20, sorted_colors(i), 'filled');
    
    % (4) 刷新图像
    % 如果暂停时间太短，我们必须手动刷新，否则看不见
    if current_pause <= min_viz_pause
        drawnow; 
        % drawnow limitrate; % (如果点太多导致卡顿，可以尝试这个)
    end
end

viz_elapsed_time = toc; % 结束计时
hold off;
disp('回放完成。');
disp(['实际回放耗时: ', num2str(viz_elapsed_time), ' 秒 (理论目标: ', num2str(target_viz_duration), ' 秒)']);




%% 其它图像绘制

signal_length = 1024;
r_loction_yld = 5.2e8+91404+1830620+703260+7049500;
ch1_yld = read_signal('..\\2023\\20230718175104.9180CH1.dat',signal_length,r_loction_yld);
bp_filtered_yld = filter_bp(ch1_yld,30e6,80e6,5);
figure
plot(bp_filtered_yld);
ylim([-80, 80]);
yticks(-80:20:80); 
xlim([0, 1024]);
xticks(0:100:1024);
xlabel('采样点');ylabel('幅值');
plot_signal_spectrum(bp_filtered_yld);




signal_length = 2e6;
% r_loction = 4.2e8;
ch1_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,3.68e8);
ch1_chj = read_signal('..\\2024 822 85933.651462CH1.dat',signal_length,3.65e8+34236156);
figure;
subplot(2,1,1);plot(ch1_yld);title('yld');
subplot(2,1,2);plot(ch1_chj);title('chj');


% 加载 .mat 文件
load('E:\双站\结果_自动校准\all_match_results_3.6-4.8_new.mat');
% 将结构体数组转换为表格
resultsTable = struct2table(all_match_results);
% 保存表格为 Excel 文件
writetable(resultsTable, 'E:\\双站\\结果_自动校准\\output.xlsx');







subplot(3,1,1);plot(filtered_signal1);title('ch1');xlabel('采样点数');ylabel('幅值');
subplot(3,1,2);plot(filtered_signal2);title('ch2');xlabel('采样点数');ylabel('幅值');
subplot(3,1,3);plot(filtered_signal3);title('ch3');xlabel('采样点数');ylabel('幅值');


plot(ch3,'b');
hold on;
plot(filtered_signal3 +50,'r');
legend('原信号','滤波后');
xlabel('采样点数');
ylabel('幅值');


[R1_x, R1_y, R1_z] = sph2cart(deg2rad(181.588691),deg2rad(49.691292),1);
% 计算仰角（单位：度）
theta = asin(R1_z) * (180/pi); % 取值范围为 [-90, 90] 度
% 计算方位角（单位：度）
phi = atan2(R1_y, R1_x) * (180/pi); % 取值范围为 [-180, 180] 度
% 如果需要将方位角转换为 [0, 360) 范围
if phi < 0
    phi = phi + 360;
end
% 输出结果
fprintf('方位角 (Azimuth): %.2f 度\n', phi);
fprintf('仰角 (Elevation): %.2f 度\n', theta);


chj_noise1 = read_signal('..\\2024 822 85933.651462CH1.dat',1e8,1e8);
chj_ch1 = read_signal('..\\2024 822 85933.651462CH1.dat',2e7,3.5e8);
data1 = load('chj_noise1.mat');
chj_Noise1 = (data1.chj_noise1)';
chj_R1 = cov(chj_Noise1); 
Q = 0.5;
%滤波
kalmanfiltered_chj1 = KalmanFilter(chj_ch1,Q,chj_R1);
subplot(2,1,1);plot(chj_ch1);title('chj');xlabel('采样点数');ylabel('幅值');
subplot(2,1,2);plot(kalmanfiltered_chj1);title('filtered_chj');xlabel('采样点数');ylabel('幅值');



signal_length = 3e7;
r_loction_yld = 3.5e8;
ch1_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction_yld);

bp_filtered_yld = filter_bp(ch1_yld,30e6,80e6,5);
plot_signal_spectrum(bp_filtered_yld);

plot_signal_spectrum(ch1_yld);
figure
plot(bp_filtered_yld);
figure
plot(ch1_yld);
% 绘制原始信号功率谱
plot_power_spectrum(ch1_yld);
% 绘制滤波后信号功率谱
plot_power_spectrum(bp_filtered_yld);


signal_length = 1.5e6;
r_loction_yld = 4.94e8;
ch1_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction_yld);
bp_filtered_yld = filter_bp(ch1_yld,30e6,80e6,5);
plot_signal_spectrum(bp_filtered_yld);

plot_signal_spectrum(ch1_yld);
figure
plot(bp_filtered_yld);
figure
plot(ch1_yld);
% 绘制原始信号功率谱
plot_power_spectrum(ch1_yld);
% 绘制滤波后信号功率谱
plot_power_spectrum(bp_filtered_yld);


signal_length = 1024;
r_loction_yld = 3.8e8+3500;
ch1_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction_yld);
bp_filtered_yld = filter_bp(ch1_yld,30e6,80e6,5);
plot_signal_spectrum(bp_filtered_yld);
figure
plot(bp_filtered_yld);


signal_length = 1024;
r_loction_yld = 3.66e8+144000;
ch1_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction_yld);
bp_filtered_yld = filter_bp(ch1_yld,30e6,80e6,5);
plot_signal_spectrum(bp_filtered_yld);
figure
plot(bp_filtered_yld);



plot_signal_spectrum(ch1_yld);

figure
plot(ch1_yld);
% 绘制原始信号功率谱
plot_power_spectrum(ch1_yld);
% 绘制滤波后信号功率谱
plot_power_spectrum(bp_filtered_yld);


signal_length = 1024;
r_loction_yld = 368631876;
ch1_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction_yld);
bp_filtered_yld = filter_bp(ch1_yld,30e6,80e6,5);
plot_signal_spectrum(bp_filtered_yld);
figure
plot(bp_filtered_yld);


figure
plot(ch1_yld);

plot_signal_spectrum(ch1_yld);


% 绘制原始信号功率谱
plot_power_spectrum(ch1_yld);
% 绘制滤波后信号功率谱
plot_power_spectrum(bp_filtered_yld);




logicalIndex =  abs(result1.t123) < 1  & abs(result1.Rcorr) > 0.65 &  result1.Start_loc < 4e8 & result1.Start_loc > 3.8e8;
filteredTable1 = result1(logicalIndex, :);
Start_loc = filteredTable1.Start_loc;
%动态图
% 归一化颜色值
Start_loc_min = min(Start_loc);
Start_loc_max = max(Start_loc);
colorValues = (Start_loc - Start_loc_min) / (Start_loc_max - Start_loc_min);

% 创建图形
figure;
hold on;
grid on;
xlabel('方位角');
xlim([0, 360]);
xticks(0:40:360);
ylabel('仰角');
ylim([0, 90]);
yticks(0:10:90);
title('目标点动态呈现');
colormap('hsv');
h_bar = colorbar;
ylabel(h_bar, '归一化起始位置');
caxis([0, 1]);

% 将颜色值划分为若干个区间（批次）
numBatches = 500;  % 控制动态速度，越小越快
bins = linspace(0, 1, numBatches + 1);

for b = 1:numBatches
    % 找到当前批次中的点
    idx = colorValues >= bins(b) & colorValues < bins(b+1);
    if any(idx)
        az = filteredTable1.Azimuth(idx);
        el = filteredTable1.Elevation(idx);
        colors = colorValues(idx);

        % 一次性绘制当前批次所有点
        scatter(az, el, 1, colors, 'filled');
        drawnow;  % 每批次刷新一次
    end
end

hold off;
disp('快速动态绘制完成。');








signal_length = 1024;
r_loction_yld = 3.8e8+3500;
ch1_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction_yld);
bp_filtered_yld = filter_bp(ch1_yld,30e6,80e6,5);
figure
plot(bp_filtered_yld);
ylim([-80, 80]);
yticks(-80:20:80); 
xlim([0, 1024]);
xticks(0:100:1024);
xlabel('采样点');ylabel('幅值');
plot_signal_spectrum(bp_filtered_yld);

% plot_signal_spectrum(bp_filtered_yld);



figure
plot(ch1_yld);
ylim([-100, 700]);
yticks(-100:100:700);
% 绘制原始信号功率谱
plot_power_spectrum(ch1_yld);
% 绘制滤波后信号功率谱
plot_power_spectrum(bp_filtered_yld);



filtered_signal1 = emd_bandpass_filter(ch1, fs, 40e6, 70e6);
figure
plot(filtered_signal1);
filtered_signal2 = emd_bandpass_filter(ch2, fs, 20e6, 90e6);
figure
plot(filtered_signal2);
filtered_signal3 = emd_bandpass_filter(ch3, fs, 20e6, 90e6);
figure
plot(filtered_signal3);



plot(segment_ch1)
ylim([-80, 80]);
yticks(-80:20:80); 
xlim([0, 1024]);
xticks(0:100:1024);
xlabel('采样点');ylabel('幅值');
figure
plot(windowed_segment_ch1)
ylim([-80, 80]);
yticks(-80:20:80); 
xlim([0, 1024]);
xticks(0:100:1024);
xlabel('采样点');ylabel('幅值');





signal_length = 0.35e8;
r_loction_yld = 3.65e8;
ch1_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction_yld);
bp_filtered_yld = filter_bp(ch1_yld,30e6,80e6,5);
figure
plot(bp_filtered_yld);
ylim([-80, 80]);
yticks(-80:20:80); 
xlim([0, 1024]);
xticks(0:100:1024);
xlabel('采样点');ylabel('幅值');
plot_signal_spectrum(bp_filtered_yld);


%% 绘制归一化后的快电场信号
signal_length = 2e4;
r_loction_yld = 380162704;
ch4_yld = read_signal('..\\20240822165932.6610CH4.dat', signal_length, r_loction_yld);

% --------------- 新增：信号归一化处理 ---------------
% 计算信号的最大值和最小值（用于归一化）
y_abs_max = max(abs(ch4_yld));  % 取信号绝对值的最大值
if y_abs_max ~= 0
    ch1_normalized = ch4_yld / y_abs_max;  % 归一化到[-1,1]范围
else
    ch1_normalized = ch4_yld;
end
baseline = movmedian(ch1_normalized, 1024);
E_fast = ch1_normalized - baseline;
% 1. 创建原始数据点的索引向量
x_indices = r_loction_yld : r_loction_yld + signal_length - 1;

% 2. 将索引转换为毫秒（ms）
time_ms = x_indices * (5 / 1e6);

% 3. 绘制归一化后的信号
% plot(time_ms, E_fast);
plot( ch1_normalized);
xlabel('时间 (ms)');
ylabel('归一化电场强度');  % 标注y轴为归一化后的幅值
title('归一化快电场信号波形');  % 可选：添加标题

%%

signal_length = 13731;
r_loction_yld = 365648853+108054;
ch4_yld = read_signal('..\\20240822165932.6610CH4.dat', signal_length, r_loction_yld);
ch1_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction_yld);
bp_filtered_ch1 = filter_bp(ch1_yld,30e6,80e6,5);
x_indices = 0 :  signal_length - 1;
time_ms = x_indices * (5 / 1e3);
baseline = movmedian(ch4_yld, 1024);
E_fast = ch4_yld - baseline;
figure
% plot(bp_filtered_ch1+300)
% hold on 
% plot(E_fast)
plot(time_ms,bp_filtered_ch1+300)
hold on 
plot(time_ms,E_fast)
title(sprintf('%d + %d 电场', r_loction_yld, signal_length));
xlabel('时间 (us)');






%%
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