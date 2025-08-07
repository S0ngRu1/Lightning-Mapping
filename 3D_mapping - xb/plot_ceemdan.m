clear; 
clc; 
close all;

%% 1. 参数定义

% ---- 信号与滤波参数 ----
fs = 200e6;  % <--- 修改这里为您正确的采样频率
freq_lower = 45e6;
freq_upper = 80e6;

% ---- CEEMDAN 函数输入参数 ----
Nstd = 0.2;         % 添加噪声的标准差 (可调)
NR = 5;           % 集合平均次数
MaxIter = 10;     % EMD内部分解最大迭代次数
SNRFlag = 1;        % 噪声模式 (1=改进版, 2=原版)，推荐使用 1

%% 2. 加载您的信号
fprintf('正在加载信号...\n');
signal_length = 1024;
r_loction_yld = 3.985e8;
try
    ch1_yld = read_signal('..\\20240822165932.6610CH2.dat', signal_length, r_loction_yld);
catch ME
    error('无法加载信号，请确保 read_signal 函数可用且文件路径正确。\n原始错误信息: %s', ME.message);
end
bp_filtered_yld = filter_bp(ch1_yld,40e6,80e6,5);
signal = bp_filtered_yld(:); % 确保是列向量
t = (0:length(signal)-1) / fs;
fprintf('信号加载成功！\n');

%% 3. ceemdan.m 进行分解
fprintf('正在执行CEEMDAN 分解...\n');
% 调用函数，并获取返回的 modes 矩阵
modes = ceemdan(signal, Nstd, NR, MaxIter, SNRFlag);

imfs_and_res = modes.'; 
fprintf('CEEMDAN 分解完成！\n');

%% 4. 筛选主导频率在指定频段内的IMF
num_total_components = size(imfs_and_res, 2);
num_imfs = num_total_components - 1; % 最后一列是残差
dominant_freqs = zeros(1, num_total_components);
selected_imfs_indices = [];
selected_imfs = [];

fprintf('正在筛选IMF...\n');
% 循环遍历所有IMF（不包括最后的残差）
for i = 1:num_imfs 
    imf_current = imfs_and_res(:, i);
    N = length(imf_current);
    
    Y = fft(imf_current);
    f_axis = (0:N-1)*(fs/N);
    half_N = ceil(N/2);
    magnitude = abs(Y(1:half_N));
    freq_axis = f_axis(1:half_N);
    
    [~, max_idx] = max(magnitude);
    dominant_freq = freq_axis(max_idx);
    dominant_freqs(i) = dominant_freq;
    
    if dominant_freq >= freq_lower && dominant_freq <= freq_upper
        selected_imfs(:, end+1) = imf_current;
        selected_imfs_indices(end+1) = i;
        fprintf('  - IMF %d (主导频率 %.2f MHz) 被选中。\n', i, dominant_freq/1e6);
    else
        fprintf('  - IMF %d (主导频率 %.2f MHz) 被舍弃。\n', i, dominant_freq/1e6);
    end
end
% 计算残差的主导频率
residual = imfs_and_res(:, end);
N_res = length(residual);
Y_res = fft(residual);
f_axis_res = (0:N_res-1)*(fs/N_res);
[~, max_idx_res] = max(abs(Y_res(1:ceil(N_res/2))));
dominant_freqs(end) = f_axis_res(max_idx_res);
fprintf('  - 残差项 (主导频率 %.2f MHz) 被舍弃。\n', dominant_freqs(end)/1e6);

%% 5. 重构信号
if isempty(selected_imfs)
    filtered_signal = zeros(size(signal));
    warning('没有找到主导频率在 %.2f MHz ~ %.2f MHz 范围内的IMF。');
else
    filtered_signal = sum(selected_imfs, 2);
end
fprintf('信号重构完成。\n');

%% 6. 结果可视化
% ---- 可视化 1: 时域信号对比 ----
figure('Name', '时域信号对比 ', 'NumberTitle', 'off');
plot(t*1e6, detrend(ch1_yld), 'k-', 'DisplayName', '原始信号');
hold on;
plot(t*1e6, filtered_signal, 'r-', 'LineWidth', 1.5, 'DisplayName', 'CEEMDAN重构信号');
hold off;
grid on;
xlabel('时间 (us)');
ylabel('幅值');
title("原始信号与CEEMDAN 带通重构信号对比");
legend('show');

% ---- 可视化 2: 频域信号对比 ----
figure('Name', '频域信号对比 ', 'NumberTitle', 'off');
% ... (此处省略与之前完全相同的绘图代码以节省篇幅，您可以直接复制过来使用)
N_fft = length(signal);
f_fft = (0:N_fft-1)*(fs/N_fft);
Y_orig = abs(fft(ch1_yld));
Y_filt = abs(fft(filtered_signal));
plot(f_fft(1:N_fft/2)/1e6, Y_orig(1:N_fft/2), 'k-', 'DisplayName', '原始信号频谱');
hold on;
plot(f_fft(1:N_fft/2)/1e6, Y_filt(1:N_fft/2), 'r-', 'LineWidth', 1.5, 'DisplayName', '重构信号频谱');
xline(freq_lower/1e6, 'b--', '频率下限');
xline(freq_upper/1e6, 'b--', '频率上限');
hold off;
grid on; xlabel('频率 (MHz)'); ylabel('幅值'); title("原始信号与重构信号的频谱 ");
legend('show'); xlim([0, fs/2/1e6]);

% ---- 可视化 3: 各IMF及残差的频谱与筛选过程 ----
figure('Name', 'IMF 频谱与筛选 ', 'NumberTitle', 'off', 'WindowState', 'maximized');
% ... (此处省略与之前完全相同的绘图代码以节省篇幅，您可以直接复制过来使用)
for i = 1:num_total_components
    subplot(ceil(num_total_components/3), 3, i);
    component = imfs_and_res(:, i); N = length(component); f = (0:N-1)*(fs/N);
    magnitude = abs(fft(component)); plot(f(1:N/2)/1e6, magnitude(1:N/2)); hold on;
    line_label = sprintf('主导: %.1f MHz', dominant_freqs(i)/1e6);
    xline(dominant_freqs(i)/1e6, 'r--', line_label, 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
    xline(freq_lower/1e6, 'k:', 'LineWidth', 1); xline(freq_upper/1e6, 'k:', 'LineWidth', 1);
    if i <= num_imfs, title_str = sprintf('IMF %d', i); else, title_str = '残差 (Residual)'; end
    if ismember(i, selected_imfs_indices), title(title_str, 'Color', 'b', 'FontWeight', 'bold');
    else, title(title_str, 'Color', 'k'); end
    xlabel('频率 (MHz)'); ylabel('幅值'); grid on; xlim([0, fs/2/1e6]);
end
sgtitle('所有IMF及残差的频谱分析 (蓝色标题表示该IMF被选中)', 'FontSize', 14, 'FontWeight', 'bold');