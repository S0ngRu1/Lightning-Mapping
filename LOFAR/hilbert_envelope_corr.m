% MATLAB code to estimate time delay between two signals based on envelope correlation

% 生成两个信号 (可以替换为实际信号)
Fs = 1000;  % 采样频率
t = 0:1/Fs:1;  % 时间向量

% 信号1：正弦信号
x1 = sin(2*pi*50*t);  

% 信号2：延迟后的正弦信号，假设延迟 0.002 秒
delay = 0.002;  
x2 = sin(2*pi*50*(t - delay));

% 添加噪声（可选）
x1 = x1 + 0.1*randn(size(x1));  
x2 = x2 + 0.1*randn(size(x2));  

% Step 1: 包络检测 (使用 Hilbert 变换)
env1 = abs(hilbert(x1));  % 信号1的包络
env2 = abs(hilbert(x2));  % 信号2的包络

% Step 2: 计算包络之间的互相关
[xcorr_env, lags] = xcorr(env1, env2);

% Step 3: 找到互相关最大值对应的时延
[~, I] = max(xcorr_env);
time_delay = lags(I) / Fs;  % 计算时间延迟

% 显示结果
fprintf('估计的时间延迟为: %f 秒\n', time_delay);

% 绘制结果
figure;
subplot(3,1,1);
plot(t, x1);
title('信号1');

subplot(3,1,2);
plot(t, x2);
title('信号2');

subplot(3,1,3);
plot(lags/Fs, xcorr_env);
title('包络互相关');
xlabel('时延 (秒)');
ylabel('互相关值');
