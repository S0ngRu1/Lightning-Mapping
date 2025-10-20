%%  静态图绘制
% --- 1. 数据准备  ---
filename = 'results\result_yld_3.5e8_4e8_window_512_128_去零飘_滤波_加窗_阈值15_30_80.txt';

% 2. 使用 readtable 函数读取数据
%    该函数会自动将第一行作为表头，并根据空格分隔各列
result1 = readtable(filename);
logicalIndex =  abs(result1.t123) < 1  & abs(result1.Rcorr) > 0.65 &  result1.Start_loc < 3.8e8 & result1.Start_loc > 3.65e8;
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
colormap('parula'); % 保持专业的颜色映射
h = colorbar;

% 颜色条标签和刻度颜色为深色
ylabel(h, '归一化发展时间', 'FontSize', 11, 'Color', 'k');
set(h, 'Color', [0.2 0.2 0.2]); % 颜色条的刻度字体颜色（深灰色）

caxis([0, 1]); % 保持颜色范围

% --- 6. 网格和整体风格 ---
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3, 'Box', 'on'); % 浅色背景下降低网格透明度



%%  动态图绘制
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
% xlim([120, 200]);
% xticks(120:20:200);
ylabel('仰角');
% ylim([10, 70]);
% yticks(10:10:70);
xlim([120, 220]);
xticks(120:20:220);
ylim([5, 85]);
yticks(5:10:85);
title('目标点动态呈现');
colormap('hsv');
h_bar = colorbar;
ylabel(h_bar, '归一化起始位置');
caxis([0, 1]);

% 将颜色值划分为若干个区间（批次）
numBatches = 3000;  % 控制动态速度，越小越快
bins = linspace(0, 1, numBatches + 1);

for b = 1:numBatches
    % 找到当前批次中的点
    idx = colorValues >= bins(b) & colorValues < bins(b+1);
    if any(idx)
        az = filteredTable1.Azimuth(idx);
        el = filteredTable1.Elevation(idx);
        colors = colorValues(idx);

        % 一次性绘制当前批次所有点
        scatter(az, el, 3, colors, 'filled');
        drawnow;  % 每批次刷新一次
    end
end

hold off;
disp('快速动态绘制完成。');

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