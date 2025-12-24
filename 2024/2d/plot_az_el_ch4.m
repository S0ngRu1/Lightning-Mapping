%% ================== 第 1 部分：加载和处理数据==================
% --- 快电场信号 ---
signal_length = 3.7e5;
r_loction_yld = 4.008e8 ;
% ch4_yld = read_signal('..\20240822165932.6610CH4.dat', signal_length, r_loction_yld);
filename = 'results\20230618125747.5480_400000000_99999999_1024_256_8_gage-20230306.txt';
result1 = readtable(filename);

point_size = 15;
sampling_interval_ns = 5;
ns_to_us = 1e-3; 
time_conversion_factor = sampling_interval_ns * ns_to_us; % us/sample

% X轴时间从 0 开始
x_indices = 0 : signal_length - 1;
time_us = x_indices * time_conversion_factor; 

% --- 1b. 准备VHF辐射源数据 ---
filter_loc_min = r_loction_yld;
filter_loc_max = r_loction_yld+signal_length;
logicalIndex =  abs(result1.t123) < 1  & ...
                abs(result1.Rcorrn) > 0.3 & ...
                result1.Start_loc < filter_loc_max & ...
                result1.Start_loc > filter_loc_min & ...
                result1.Elevation < 80 & ...
                result1.Azimuth > 255 & ...
                result1.Azimuth < 330 ;
            
filteredTable1 = result1(logicalIndex, :);

% 将VHF样本位置转换为相对时间 (us)，使其从 0 开始
vhf_time_us = (filteredTable1.Start_loc - r_loction_yld) * time_conversion_factor;
vhf_elevation = filteredTable1.Elevation;

figure('Color', [1 1 1], 'Position', [100, 100, 900, 600]);
ax = gca; % 获取当前坐标轴

% 计算X轴范围 
time_min_us = 0;
time_max_us = (signal_length - 1) * time_conversion_factor;

% 直接绘制散点图
scatter(ax, vhf_time_us, vhf_elevation, ...
        point_size, ...
        vhf_time_us, ... % 颜色映射数据
        'filled', 'MarkerFaceAlpha', 0.7);

ylabel(ax, '仰角 (Elevation / °)', 'FontSize', 12, 'FontWeight', 'bold');
title('VHF辐射源仰角时序图', 'FontSize', 16, 'FontWeight', 'bold');

% 设置Y轴颜色为黑色
set(ax, 'YColor', [0 0 0]);

% X轴设置
xlabel(ax, '时间 (us)', 'FontSize', 12, 'FontWeight', 'bold');
xticks(0:100:time_max_us);

% 颜色条设置 
h = colorbar(ax);
ylabel(h, '时间 (us)', 'FontSize', 11);
colormap(ax, 'jet'); 
clim(ax, [time_min_us, time_max_us]); 

% 网格和框线
grid(ax, 'on');
set(ax, 'FontSize', 11, 'LineWidth', 1.2, 'Box', 'on', 'GridAlpha', 0.3, 'GridLineStyle', '--');