% ========================================================================
% MATLAB 代码：证明 T 在不同信噪比下的边际收益差异
% ========================================================================

clear; close all; clc;

%% 1. 参数设置
f1 = 1e6;
f2 = 5e6;
Constant_Term = 3 / (8 * pi^2 * (f2^3 - f1^3));

% 定义时间轴：0 到 25000 ns
t_vec = linspace(100e-9, 25000e-9, 1000); 

% 定义两个典型的信噪比场景
SNR_low = 1.0;   % 低信噪比 (e.g., 0 dB 左右)
SNR_high = 10.0; % 高信噪比 (e.g., 20 dB 左右)

%% 2. 计算方差
% 公式: sigma^2 ~ (1+2SNR)/(SNR^2 * T)
sigma_low = Constant_Term * (1 + 2*SNR_low) ./ (SNR_low^2 * t_vec);
sigma_high = Constant_Term * (1 + 2*SNR_high) ./ (SNR_high^2 * t_vec);

%% 3. 绘图 (核心部分：使用双纵轴展示对比)
figure('Color', 'w', 'Position', [100, 100, 800, 500]);

% --- 左纵轴：低信噪比情况 ---
yyaxis left
plot(t_vec*1e9, sigma_low, 'b-', 'LineWidth', 2.5);
ylabel(['低信噪比 (SNR=' num2str(SNR_low) ') 的方差 \sigma^2'], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'YColor', 'b'); % 设置左轴颜色为蓝色
% 标记显著下降区
text(2000, max(sigma_low)*0.6, {'\leftarrow 增加 T 带来', '显著精度提升'}, 'Color', 'b', 'FontSize', 11);

hold on;

% --- 右纵轴：高信噪比情况 ---
yyaxis right
plot(t_vec*1e9, sigma_high, 'r--', 'LineWidth', 2.5);
ylabel(['高信噪比 (SNR=' num2str(SNR_high) ') 的方差 \sigma^2'], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'YColor', 'r'); % 设置右轴颜色为红色
ylim([0, max(sigma_high)*1.2]); % 调整右轴范围，避免曲线贴底

% --- 添加辅助线和标注 ---
grid on;
xlabel('积分时间 T (ns)', 'FontSize', 12, 'FontWeight', 'bold');
title('论证：不同信噪比下增加积分时间 T 的收益对比', 'FontSize', 14);

% 绘制“饱和/无效区”分界线 (例如在 10000ns 处)
xline(10000, 'k:', 'LineWidth', 2);
text(10500, max(sigma_high)*0.8, {'收益饱和区', '继续增加 T 效果不明显'}, 'FontSize', 11, 'BackgroundColor', 'w');

% 添加图例
legend({'低 SNR 曲线 (左轴)', '高 SNR 曲线 (右轴)'}, 'Location', 'northeast', 'FontSize', 11);

%% 4. 导数（边际收益）量化输出
% 为了让你的证明更有力，计算一下在 T=1000ns 处的下降速率
rate_low = abs(diff(sigma_low(1:2)) / diff(t_vec(1:2)));
rate_high = abs(diff(sigma_high(1:2)) / diff(t_vec(1:2)));

fprintf('=== 证明数据 (用于论文描述) ===\n');
fprintf('在 T 较小时 (起始阶段)：\n');
fprintf('低 SNR 下，每增加单位时间，方差降低速率(收益)为: %.2e\n', rate_low);
fprintf('高 SNR 下，每增加单位时间，方差降低速率(收益)为: %.2e\n', rate_high);
fprintf('结论：低 SNR 下的改进速度是高 SNR 下的 %.1f 倍。\n', rate_low/rate_high);