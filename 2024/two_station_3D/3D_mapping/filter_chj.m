
% 使用实际的闪电信号数据
signal_length = 6000;
start_read_loc_chj = 4.5e8;
chj3 = read_signal('../2024 822 85933.651462CH3.dat',signal_length,start_read_loc_chj+165/5);

% 创建时间向量
fs = 200e6;  
t = (0:length(chj3)-1)/fs;

noisy_signal = normalize(chj3);  % 含噪声原始信号
%% 不同滤波器尝试

% % 1. 带通滤波
filtered_signal1 = normalize(filter_bp(chj3, 20e6 ,80e6 ,5));

% 2. 中值滤波
window_size = 5;  % 可调整窗口大小
median_filtered = normalize(medfilt1(noisy_signal, window_size));

% 3. 小波去噪
level = 4;  % 分解层数
wavelet = 'db4';  % 小波基
[c,l] = wavedec(noisy_signal, level, wavelet);
% 自适应阈值去噪
sigma = median(abs(c))/0.6745;  % 估计噪声标准差
thr = sigma*sqrt(2*log(length(noisy_signal)));
c_denoised = wthresh(c, 's', thr);
wavelet_filtered = normalize(waverec(c_denoised, l, wavelet));


% 4. 零相位滤波器
corner_frequency = 100e6; % 角频率 (Hz)
% 应用零相位滤波器
zero_phase_filtered = normalize(zero_phase_filter(noisy_signal, fs, corner_frequency));


% 5. 卡尔曼滤波器
% 双指数脉冲参数
A0 = 6;                   % 最大幅值
tau1 = 5;                % 主衰减时间常数500μs
tau2 = 5;                % 快速上升时间常数10μs（陡峭上升）
% 初始状态估计（可基于先验知识调整）
x0 = [A0; tau1; tau2];      % [初始幅度; τ1; τ2]
% 过程噪声协方差（控制模型信任度）
Q = diag([1e3, 1, 1]); % [幅度噪声, τ1噪声, τ2噪声]
% 观测噪声协方差（根据实际噪声方差设置）
R = 1000;            % R = 噪声方差 = 10^2 = 100
% 初始估计协方差（反映初始状态不确定性）
P0 = diag([1e4, 1e4, 1e4]); 
% 执行滤波
[ekf_filtered, x_est] = ekf_lightning_pulse(noisy_signal, fs, x0, Q, R, P0);
ekf_filtered = normalize(ekf_filtered); % 归一化滤波结果


% 绘图
figure('Position', [100 100 1200 800]);

% 带通滤波结果
subplot(4,2,1);
plot(t*1e6, noisy_signal, 'b', 'LineWidth', 0.5);
hold on;
plot(t*1e6, filtered_signal1, 'r', 'LineWidth', 1);
title('带通滤波结果');
xlabel('时间 (μs)');
ylabel('幅度');
legend('原始信号', '滤波后');
grid on;

% 中值滤波结果
subplot(4,2,2);
plot(t*1e6, noisy_signal, 'b', 'LineWidth', 0.5);
hold on;
plot(t*1e6, median_filtered, 'r', 'LineWidth', 0.5);
title('中值滤波结果');
xlabel('时间 (μs)');
ylabel('幅度');
legend('原始信号', '滤波后');
grid on;

% 小波去噪结果
subplot(4,2,3);
plot(t*1e6, noisy_signal, 'b', 'LineWidth', 0.5);
hold on;
plot(t*1e6, wavelet_filtered, 'r', 'LineWidth',  0.5);
title('小波去噪结果');
xlabel('时间 (μs)');
ylabel('幅度');
legend('原始信号', '滤波后');
grid on;

% 零相位滤波结果
subplot(4,2,4);
plot(t*1e6, noisy_signal, 'b', 'LineWidth', 0.5);
hold on;
plot(t*1e6, zero_phase_filtered, 'r', 'LineWidth',  0.5);
title('零相位滤波结果');
xlabel('时间 (μs)');
ylabel('幅度');
legend('原始信号', '滤波后');
grid on;

% EKF结果
subplot(4,2,5);
plot(t*1e6, noisy_signal, 'b', 'LineWidth', 0.5);
hold on;
plot(t*1e6, ekf_filtered, 'r', 'LineWidth',  0.5);
title('卡尔曼滤波结果');
xlabel('时间 (μs)');
ylabel('幅度');
legend('原始信号', '滤波后');
grid on;

% 原信号
subplot(4,2,6);
plot(t*1e6, noisy_signal, 'b', 'LineWidth', 0.5);
title('原信号');
xlabel('时间 (μs)');
ylabel('幅度');
grid on;
% 所有滤波结果对比
subplot(4,2,[7,8]);
plot(t*1e6, noisy_signal, 'k', 'LineWidth', 0.5);
hold on;
plot(t*1e6, filtered_signal1, 'r', 'LineWidth',  0.5);
plot(t*1e6, median_filtered, 'g', 'LineWidth',  0.5);
plot(t*1e6, wavelet_filtered, 'b', 'LineWidth',  0.5);
plot(t*1e6, zero_phase_filtered, 'r', 'LineWidth',  0.5);
title('所有方法对比');
xlabel('时间 (μs)');
ylabel('幅度');
legend('原始信号', '带通滤波', '中值滤波', '小波去噪');
grid on;


%% 小波滤波的改进
% 基本参数设置
% 原始信号

% 基本参数设置
level = 4;          % 分解层数
wavelet = 'db4';    % 小波基 

% 1. 改进的小波去噪 - SURE阈值
%基于Stein无偏风险估计原理的自适应阈值方法。它通过对给定阈值的似然估计，再将非似然T最小化，从而得到所选的阈值。
[c, l] = wavedec(noisy_signal, level, wavelet);
thr = zeros(level, 1);
c_denoised = c;
len = length(noisy_signal);

% 分层处理系数
last = 0;
for i = 1:level
    first = last + 1;
    last = first + l(i) - 1;
    
    % 每层系数
    d = c(first:last);
    
    % SURE阈值
    N = length(d);
    sigma = median(abs(d))/0.6745;
    thr(i) = sigma*sqrt(2*log(N));
    
    % 改进的软阈值函数
    idx = abs(d) <= thr(i);
    d(idx) = 0;
    d(~idx) = sign(d(~idx)).*(abs(d(~idx)) - thr(i));
    
    c_denoised(first:last) = d;
end

% 重构信号
denoised_signal1 = normalize(waverec(c_denoised, l, wavelet));

% 2. 改进的小波去噪 - 自适应阈值
%自适应阈值方法根据信号的特点和噪声的特点来选择适当的阈值，并通过构建适应性的阈值函数进行阈值的确定
[c2, l2] = wavedec(noisy_signal, level, wavelet);
c_denoised2 = c2;
last = 0;

for i = 1:level
    first = last + 1;
    last = first + l2(i) - 1;
    d = c2(first:last);
    
    % 自适应阈值计算
    sigma = median(abs(d))/0.6745;
    N = length(d);
    thr_univ = sigma*sqrt(2*log(N));  % 通用阈值
    thr_mini = sigma*sqrt(2*log(log(N))); % minimax阈值
    
    % 混合阈值策略
    energy = sum(d.^2);
    if energy > median(d.^2)*N
        % 信号能量大，使用较小阈值
        thr_adapt = thr_mini;
    else
        % 信号能量小，使用较大阈值
        thr_adapt = thr_univ;
    end
    
    % 改进的阈值函数（平滑过渡）
    beta = 1.5;  % 平滑因子
    d_processed = sign(d).*max(0, abs(d) - thr_adapt./(1 + exp(-beta*(abs(d)-thr_adapt))));
    
    c_denoised2(first:last) = d_processed;
end

% 重构信号
denoised_signal2 = normalize(waverec(c_denoised2, l2, wavelet));

% 3. 相位保持小波去噪
%在去噪过程中保持信号的相位信息
[c3, l3] = wavedec(noisy_signal, level, wavelet);
c_denoised3 = c3;
last = 0;

for i = 1:level
    first = last + 1;
    last = first + l3(i) - 1;
    d = c3(first:last);
    
    % 计算局部方差
    wind_size = min(length(d), 32);
    local_var = zeros(size(d));
    for j = 1:length(d)
        wind_start = max(1, j - floor(wind_size/2));
        wind_end = min(length(d), j + floor(wind_size/2));
        local_var(j) = var(d(wind_start:wind_end));
    end
    
    % 自适应阈值
    sigma = median(abs(d))/0.6745;
    thr_base = sigma*sqrt(2*log(length(d)));
    
    % 基于局部方差的阈值调整
    thr_local = thr_base ./ (1 + local_var/max(local_var));
    
    % 软阈值处理
    for j = 1:length(d)
        if abs(d(j)) <= thr_local(j)
            d(j) = 0;
        else
            d(j) = sign(d(j))*(abs(d(j)) - thr_local(j));
        end
    end
    
    c_denoised3(first:last) = d;
end

% 重构信号
denoised_signal3 = normalize(waverec(c_denoised3, l3, wavelet));

% 4. 标准小波去噪
%使用固定的阈值或基于噪声估计的阈值进行去噪处理
[c4, l4] = wavedec(noisy_signal, level, wavelet);
sigma = median(abs(c4))/0.6745;
% thr = sigma*sqrt(2*log(length(noisy_signal)));  % 通用阈值
thr = sigma*sqrt(log(length(noisy_signal))); 
c_denoised4 = wthresh(c4, 's', thr);  % 使用小波软阈值
denoised_signal4 = normalize(waverec(c_denoised4, l4, wavelet));

% 绘图比较结果
figure('Position', [100 100 1200 800]);

% 原始信号
subplot(5,1,1);
plot(t*1e6, noisy_signal, 'b', 'LineWidth', 0.5);
title('原始闪电信号');
xlabel('时间 (μs)');
ylabel('幅度');
grid on;

% SURE阈值结果
subplot(5,1,2);
plot(t*1e6, denoised_signal1, 'b', 'LineWidth', 0.5);
title('SURE阈值去噪结果');
xlabel('时间 (μs)');
ylabel('幅度');
grid on;

% 自适应阈值结果
subplot(5,1,3);
plot(t*1e6, denoised_signal2, 'b', 'LineWidth', 0.5);
title('自适应阈值去噪结果');
xlabel('时间 (μs)');
ylabel('幅度');
grid on;

% 相位保持结果
subplot(5,1,4);
plot(t*1e6, denoised_signal3, 'b', 'LineWidth', 0.5);
title('相位保持去噪结果');
xlabel('时间 (μs)');
ylabel('幅度');
grid on;

% 标准小波去噪结果
subplot(5,1,5);
plot(t*1e6, denoised_signal4, 'b', 'LineWidth', 0.5);
title('标准小波去噪结果');
xlabel('时间 (μs)');
ylabel('幅度');
grid on;

