% 读取数据
flash_data = readtable('result_4-6e8_w1024_256_R0.8.txt', 'Whitespace', ' ');
logicalIndex = abs(flash_data.t123) < 1 & abs(flash_data.Rcorr) > 0.5;
% 假设你的数据存储在变量 data 中
filteredTable = flash_data(logicalIndex, :);
azimuthData = filteredTable.Azimuth;
elevationData = filteredTable.Elevation;

% 设置移动平均窗口大小
windowSize = 5;
smoothedAzimuth = movmean(azimuthData, windowSize);
smoothedElevation = movmean(elevationData, windowSize);

% 将平滑后的数据更新到原表格中
filteredTable.SmoothedAzimuth = smoothedAzimuth;
filteredTable.SmoothedElevation = smoothedElevation;

Start_loc = filteredTable.Start_loc;
colorValues = (Start_loc - 3e8) / 2e8;
figure;
scatter(filteredTable.SmoothedAzimuth,filteredTable.SmoothedElevation, 1, colorValues, 'filled');
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