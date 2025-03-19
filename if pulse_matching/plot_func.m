%结果1
% logicalIndex = abs(result1.t123) < 1 & abs(result1.Rcorr) > 0.3;
logicalIndex = abs(result1.t123) <1 & abs(result1.Rcorr) > 0.3 &  result1.Start_loc < 523418116 & result1.Start_loc > 451501808;
filteredTable1 = result1(logicalIndex, :);
Start_loc = filteredTable1.Start_loc;
colorValues = (Start_loc - 3e8) / 2e8;
figure;
scatter(filteredTable1.Azimuth,filteredTable1.Elevation, 1, colorValues, 'filled');
title('Azimuth vs Elevation');
xlabel('Azimuth');
xlim([0, 360]);
xticks(0:40:360);
ylabel('Elevation');
ylim([-40, 100]);
yticks(-40:20:100);
colormap('hot');
colorbar;
caxis([0, 1.5]);
grid on;

%结果2
figure;
% logicalIndex = abs(result2.t123) < 1 & abs(result2.Rcorr) > 0.3 &  result2.Start_loc > 500000000 ;
logicalIndex = abs(result2.t123) < 1 & abs(result2.Rcorr) > 0.5;
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
colormap('hsv');
colorbar;
caxis([0, 1.5]);
grid on;

%结果3
logicalIndex = abs(result3.t123) > 0.001 & abs(result3.Rcorr) < 0.2 &  result3.Start_loc < 600000000 & result3.Start_loc > 400000000;
filteredTable3 = result3(logicalIndex, :);
Start_loc = filteredTable3.Start_loc;
colorValues = (Start_loc - 3e8) / 2e8;
figure;
scatter(filteredTable3.Azimuth,filteredTable3.Elevation, 1, colorValues, 'filled');
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































% logicalIndex = abs(Untitled.VarName11) < 1 & abs(Untitled.VarName12) > 0.25 & Untitled.Paramers >500000000;
logicalIndex = abs(Untitled.VarName11) < 1 & abs(Untitled.VarName12) > 0.3 & Untitled.Paramers > 400000000 & Untitled.Paramers < 600000000 ;
filteredTable3 = Untitled(logicalIndex, :);
Start_loc = filteredTable3.Paramers;
colorValues = (Start_loc - 3e8) / 2e8;
figure;
scatter(filteredTable3.VarName5,filteredTable3.VarName6, 1, colorValues, 'filled');
title('Azimuth vs Elevation');
xlabel('Azimuth');
xlim([0, 360]);
xticks(0:40:360);
ylabel('Elevation');
ylim([-40, 100]);
yticks(-40:20:100);
colormap('hot');
colorbar;
caxis([0, 1.5]);
grid on;



% logicalIndex = abs(result1.t123) < 1 & abs(result1.Rcorr) > 0.3;
% filteredTable1 = result1(logicalIndex, :);
% index = 1:256:size(filteredTable1, 1);
% Start_loc = filteredTable1.Start_loc;
% colorValues = (Start_loc - 3e8) / 2e8;
% figure;
% scatter(filteredTable1.Azimuth,filteredTable1.Elevation, 1, colorValues, 'filled');
% title('Azimuth vs Elevation');
% xlabel('Azimuth');
% xlim([0, 360]);
% xticks(0:40:360);
% ylabel('Elevation');
% ylim([-40, 100]);
% yticks(-40:20:100);
% colormap('hot');
% colorbar;
% caxis([0, 1.5]);
% grid on;






% 绘制两个波形
plot(processed_yld_signal-50, 'b');
hold on;
plot(processed_chj_signal, 'r');
legend('yld', 'chj');
xlabel('采样点数');
ylabel('幅值');
title('bp filter');






yld_signal = read_signal('../20240822165932.6610CH1.dat',6000,start_read_loc_yld);
chj_signal = read_signal('../2024 822 85933.651462CH1.dat',6000,start_read_loc_chj);
x = downsample(yld_signal,1);
y = downsample(chj_signal,1);
subplot(2,1,1);plot(x);title('yld');xlabel('采样点数');ylabel('幅值');
subplot(2,1,2);plot(y);title('chj');xlabel('采样点数');ylabel('幅值');




scatter(filteredTable1.cos_alpha_opt,filteredTable1.cos_beta_opt,2);
xlim([-1, 1]);
% 指定y轴范围和刻度标记
ylim([-1, 1]);

scatter(Untitled.VarName5,Untitled.VarName6,1);
xlim([0, 360]);
xticks(0:40:360);
% 指定y轴范围和刻度标记
ylim([0, 90]);
yticks(0:20:90);
scatter(result1.VarName10,result1.VarName11);
scatter(result1.cos,result1.cos1);
histogram(result.VarName11);
histogram(Untitled.VarName6);
subplot(3,1,1);plot(signal3);title('ch1');xlabel('采样点数');ylabel('幅值');
%设置横坐标间隔为32

subplot(3,1,2);plot(filtered_signal3);title('ch2');xlabel('采样点数');ylabel('幅值');

subplot(3,1,3);plot(filtered_signal3);title('ch3');xlabel('采样点数');ylabel('幅值');


plot(lag12_msw,R12_msw);
plot(ch1);
 %绘制上采样对比图
plot(ch3_gcc_new, 'b');
hold on;
plot(ch3_upsp(:,1), ch3_upsp(:,2), 'r--');
legend('ch3波形', '上采样ch3波形');
xlabel('采样点数');
ylabel('幅值');
plot(lag12_gcc,R12_gcc);

 %绘制上采样对比图
plot(filtered_signal1, 'b');
hold on;
plot(ch1, 'r--');
legend('ch2波形', '上采样ch2波形');
xlabel('采样点数');
ylabel('幅值');
axis auto


% 绘制两个波形
plot(ch1, 'b');
hold on;
plot(IMF(:,3), 'r--');
legend('ch1波形', '过滤后的ch1波形');
xlabel('采样点数');
ylabel('幅值');
title('根据互相关曲线的最大值进行数据平移的结果');

% 绘制两个波形
plot(ch1_gcc_new, 'b');
hold on;
plot(ch3_gcc_new, 'r--');
legend('原信号', '平移后');
xlabel('采样点数');
ylabel('幅值');
title('结果');



num_bins = 100; % 定义直方图的柱数
[counts, edges] = histcounts(matched_max_r_gccs, num_bins);

% 计算每个区间的中心点，用于绘制频率分布图
centers = edges(1:end-1) + diff(edges) / 2;

% 绘制频率分布图
figure;
bar(centers, counts);
xlabel('相关系数大小');
ylabel('数量');
title('匹配到的相关系数分布图');
grid on; % 添加网格