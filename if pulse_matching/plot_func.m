%结果1
logicalIndex = abs(result1.t123) <1 & abs(result1.Rcorr) > 0.3 &  result1.Start_loc < 4.8e8 & result1.Start_loc > 3.6e8;
filteredTable1 = result1(logicalIndex, :);
Start_loc = filteredTable1.Start_loc;
% colorValues = (Start_loc - 3e8) / 2e8;
colorValues = (Start_loc - min(Start_loc)) / (max(Start_loc) - min(Start_loc));
% 转换 Azimuth 范围从 0-360 到 -180-180
% filteredTable1.Azimuth = mod(filteredTable1.Azimuth - 180, 360) - 180;
figure;
scatter(filteredTable1.Azimuth,filteredTable1.Elevation, 2, colorValues, 'filled');
title('Azimuth vs Elevation');
xlabel('Azimuth');
xlim([0, 360]);
xticks(0:40:360);
ylabel('Elevation');
ylim([-40, 100]);
yticks(-40:20:100);
colormap('hsv');
colorbar;
caxis([0, 1.5]);
grid on;
% 
% %结果2
% figure;
%  logicalIndex = abs(result1.t123) <1 & abs(result1.Rcorr) > 0.4 &  result1.Start_loc < 4.5e8 & result1.Start_loc > 4e8;
logicalIndex = abs(result2.t123) < 1 & abs(result2.Rcorr) > 0.3;
filteredTable2 = result2(logicalIndex, :);
index = 1:256:size(filteredTable2, 1);
Start_loc = filteredTable2.Start_loc;
colorValues = (Start_loc - 3e8) / 2e8;
figure;
scatter(filteredTable2.Azimuth,filteredTable2.Elevation, 1, colorValues, 'filled');
title('Azimuth vs Elevation');
xlabel('Azimuth');
xlim([0, 360]);
xticks(0:40:360);
ylabel('Elevation');
ylim([-40, 100]);
yticks(-40:20:100);
colormap('jet');
colorbar;
caxis([0, 1.5]);
grid on;

% %结果3
% logicalIndex = abs(result3.t123) > 0.001 & abs(result3.Rcorr) < 0.2 &  result3.Start_loc < 600000000 & result3.Start_loc > 400000000;
% filteredTable3 = result3(logicalIndex, :);
% Start_loc = filteredTable3.Start_loc;
% colorValues = (Start_loc - 3e8) / 2e8;
% figure;
% scatter(filteredTable3.Azimuth,filteredTable3.Elevation, 1, colorValues, 'filled');
% title('Azimuth vs Elevation');
% xlabel('Azimuth');
% xlim([0, 360]);
% xticks(0:40:360);
% ylabel('Elevation');
% ylim([-40, 100]);
% yticks(-40:20:100);
% colormap('hsv');
% colorbar;
% caxis([0, 1.5]);
% grid on;



data_original1=window_plus(N,data_original1);
data_original2=window_plus(N,data_original2);
Wc=2*80e6/200e6;
Wp=2*30e6/200e6;
data_filter1 = datafilter(data_original1,Wc,Wp);
data_filter2 = datafilter(data_original2,Wc,Wp);
data_noisefilter1 = datafilter(data_noise1,Wc,Wp);
data_noisefilter2 = datafilter(data_noise2,Wc,Wp);


Wc=2*80e6/200e6;
Wp=2*30e6/200e6;
N = 2e7;
r_loction = 256;
ch1 = read_signal('..\\20240822165932.6610CH1.dat',N,r_loction);
data_original1=window_plus(N,ch1);
data_filter1 = datafilter(data_original1,Wc,Wp);


signal_length = 1024;
r_loction_yld = 469401640;
r_loction_chj = 504138238;
ch_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction_yld-256);
ch_chj = read_signal('..\\2024 822 85933.651462CH1.dat',signal_length,r_loction_chj);
% filtered_yld = rfi_filter(ch_yld,1536);
% filtered_chj = rfi_filter(ch_chj,1536);
filtered_yld = filter_bp(ch_yld,20e6,80e6,5);
filtered_chj = filter_bp(ch_chj,20e6,80e6,5);
subplot(2,1,1);plot(filtered_yld);title('yld');xlabel('采样点数');ylabel('幅值');
subplot(2,1,2);plot(filtered_chj);title('chj');xlabel('采样点数');ylabel('幅值');

 filtered_chj_signal1 = waveletDenoiseAdaptive(ch2, level, wavelet);
 filtered_chj_signal2 = filter_bp(ch2,20e6,80e6,5);
x1 = downsample(ch1,50);
x2 = downsample(ch2,50);
x3 = downsample(ch3,50);
% x3 = downsample(ch3,50);
subplot(3,1,1);plot(ch1);title('ch1');xlabel('采样点数');ylabel('幅值');
subplot(3,1,2);plot(data_filter1);title('wf');xlabel('采样点数');ylabel('幅值');
subplot(3,1,3);plot(filter2);title('bp');xlabel('采样点数');ylabel('幅值');

%滑动25%滤波


%信号波形与阈值对比图
N = 5e7; % 信号总长度
subsignal_length = 4000;
M = ceil(N / subsignal_length); % 总段数
subsignal_start = 1:subsignal_length:length(filtered_signal1);
for subi = 1:numel(subsignal_start)
subsignal1 = filtered_signal1(subsignal_start(subi):subsignal_start(subi)+subsignal_length-1);
threshold =  mean(abs(subsignal1)) + 3*std(subsignal1);
end
% 假设每段的阈值已经计算好，存储在 thresholds 中
threshold = rand(1, M) * 5; % 示例阈值，实际中用你的阈值替换
% 绘制信号
figure;
plot(filtered_signal1, 'b', 'LineWidth', 1.5); % 绘制信号，蓝色线
hold on;
% 绘制每段的阈值
for i = 1:M
    start_index = (i - 1) * subsignal_length + 1;
    end_index = min(i * subsignal_length, N); % 防止最后一段超出信号长度
    x = [start_index, end_index]; % 当前段的 x 范围
    y = [threshold(i), threshold(i)]; % 当前段的阈值
    plot(x, y, 'r--', 'LineWidth', 1.5); % 用红色虚线绘制阈值
end
% 添加图例和标签
legend('filtered_signal1', 'Thresholds');
xlabel('Sample Index');
ylabel('Amplitude');
title('Signal and Thresholds Comparison');
grid on; % 添加网格
hold off;





%3D转2D
x = filtered_S(:, 1);
y = filtered_S(:, 2);
z = filtered_S(:, 3);
% 将直角坐标转换为球坐标
[azimuth, elevation, ~] = cart2sph(x, y, z);
% 将弧度转换为角度
azimuth = azimuth * 180 / pi;
elevation = elevation * 180 / pi;
% 调整方位角范围为 0-360 度
azimuth = mod(azimuth, 360);
azimuth(azimuth < 0) = azimuth(azimuth < 0) + 360;
% 调整仰角范围为 0-90 度
elevation = abs(elevation); % 确保仰角为非负值
elevation(elevation > 90) = 90; % 限制最大值为 90 度
% 绘制方位角和仰角的散点图
figure;
scatter(azimuth, elevation, 1, 'filled');
xlabel('方位角 (度)');
ylabel('仰角 (度)');
title('方位角和仰角分布');
grid on;
axis([0 360 0 90]); % 设置坐标轴范围
colorbar; % 添加颜色条（可选）




subplot(3,1,1);plot(ch1);title('ch1');xlabel('采样点数');ylabel('幅值');
subplot(3,1,2);plot(ch2);title('ch2');xlabel('采样点数');ylabel('幅值');
subplot(3,1,3);plot(ch3);title('ch3');xlabel('采样点数');ylabel('幅值');







fs = 200e6;
signal_length = 71916308;
r_loction = 451501808;
ch1_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction);
ch1_chj = read_signal('..\\2024 822 85933.651462CH1.dat',signal_length,r_loction + 3.3e7);
ch2 = read_signal('..\\2024 822 85933.651462CH2.dat',signal_length,r_loction);
ch3 = read_signal('..\\2024 822 85933.651462CH3.dat',signal_length,r_loction+165/5);

% 设计带通滤波器
[b, a] = butter(4, [20e6, 80e6] / (fs / 2), 'bandpass'); 
% 对信号进行带通滤波
filtered_signal1 = filter(b, a, ch1);
filtered_signal2 = filter(b, a, ch2);
filtered_signal3 = filter(b, a, ch3);

%bp
filtered_signal1_yld = filter_bp(ch1_yld, 20e6 ,80e6 ,5);
filtered_signal1_chj = filter_bp(ch1_chj, 20e6 ,80e6 ,5);
filtered_signal3 = filter_bp(ch3, 20e6 ,80e6 ,5);








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



