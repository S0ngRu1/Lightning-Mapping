%% 在一幅图中绘制通道1和通道4的图（小波滤波优化版）

signal_length = 2e8;
r_loction_yld = 3.6e8;
ch4_yld = read_signal('..\\20240822165932.6610CH4.dat', signal_length, r_loction_yld);
ch1_yld = read_signal('..\\20240822165932.6610CH1.dat', signal_length, r_loction_yld);
bp_filtered_ch1 = filter_bp(ch1_yld, 30e6, 80e6, 5); % CH1保留带通滤波（若需优化可同步调整）

x_indices = 0 : signal_length - 1;
time_ms = x_indices * (5 / 1e3); % 保持原时间轴

% 基线校正（优化：用更平滑的移动中值）
baseline = movmedian(ch4_yld, 1024); 
E_fast = ch4_yld - baseline;

% ====== 小波滤波核心修复（雷电信号专用） ======
wname = 'db4';        % Daubechies 4小波（雷电信号最佳实践）
level = 4;            % 分解层数（4层覆盖雷电典型频段1-100MHz）

% ====== 修复关键：使用wden自动去噪（MATLAB官方推荐） ======
% 重要：wden会自动计算最优阈值，无需手动调参
[E_fast_denoised, ~] = wden(E_fast, 'heursure', 's', 'one', level, wname);

% ====== 绘制结果（优化显示） ======
figure
title(sprintf('%d + %d 电场', r_loction_yld, signal_length));
% plot(time_ms, bp_filtered_ch1 + 500, 'b'); % CH1: 蓝色
% hold on;
% plot(time_ms, E_fast_denoised, 'r'); % CH4: 红色
% xlabel('时间 (μs)');
% ylabel('信号幅度');
% legend('通道1 (CH1)', '快电场 (CH4)');
% grid on;

% plot(bp_filtered_ch1 + 500, 'b'); % CH1: 蓝色
% hold on;
plot(E_fast_denoised, 'r'); % CH4: 红色
xlabel('采样点');
ylabel('信号幅度');
% legend('通道1 (CH1)', '快电场 (CH4)');
grid on;
