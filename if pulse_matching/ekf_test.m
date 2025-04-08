%% 闪电脉冲信号EKF滤波完整示例
% 作者：AI助手
% 功能：针对200MHz采样、最大幅值660的闪电脉冲信号进行EKF滤波

clc; clear; close all;

%% 1. 生成模拟闪电脉冲信号
fs = 200e6;                 % 采样频率200MHz


signal_length = 5E6; 
start_read_loc_yld = 451518508;
yld_signal = read_signal('../20240822165932.6610CH1.dat',signal_length,start_read_loc_yld);
n_points = signal_length;             % 总采样点数
t = (0:n_points-1)' / fs;   % 时间向量

% 双指数脉冲参数
A0 = 660;                   % 最大幅值
tau1 = 5e-4;                % 主衰减时间常数500μs
tau2 = 1e-5;                % 快速上升时间常数10μs（陡峭上升）

%% 3. 设置EKF参数并执行滤波
% 初始状态估计（可基于先验知识调整）
x0 = [A0; tau1; tau2];      % [初始幅度; τ1; τ2]

% 过程噪声协方差（控制模型信任度）
Q = diag([1e3, 1e-8, 1e-8]); % [幅度噪声, τ1噪声, τ2噪声]

% 观测噪声协方差（根据实际噪声方差设置）
R = 100;            % R = 噪声方差 = 10^2 = 100

% 初始估计协方差（反映初始状态不确定性）
P0 = diag([1e4, 1e-4, 1e-4]); 

% 执行滤波
[filtered, x_est] = ekf_lightning_pulse(yld_signal, fs, x0, Q, R, P0);

%% 4. 结果可视化
% 原始信号 vs 滤波结果
plot(t, yld_signal, 'Color', [0.7,0.7,0.7], 'LineWidth', 0.5); hold on;
plot(t, filtered-50, 'r', 'LineWidth', 0.5);
legend('含噪信号', '真实信号', 'EKF滤波结果', 'Location', 'northeast');
title('信号滤波效果对比');
xlabel('时间 (s)'); ylabel('幅度');
xlim([0, 1e-3]); % 显示前1ms细节
grid on;



%% 2. 定义EKF滤波函数
function [filtered_signal, x_estimates] = ekf_lightning_pulse(raw_signal, fs, x0, Q, R, P0)
    % EKF_LIGHTNING_PULSE 扩展卡尔曼滤波处理闪电脉冲信号（双指数模型）
    % 状态变量: x = [幅度A; 衰减时间常数τ1; 上升时间常数τ2]
    
    dt = 1/fs;              % 采样间隔
    n_samples = length(raw_signal);
    filtered_signal = zeros(n_samples, 1);
    x_estimates = zeros(3, n_samples);
    
    % 初始化
    x_hat = x0(:);          % 转为列向量
    P = P0;
    
    for k = 1:n_samples
        % ---------- 预测步骤 ----------
        A_prev = x_hat(1);
        tau1_prev = x_hat(2);
        tau2_prev = x_hat(3);
        
        % 状态预测（双指数模型）
        x_pred = [
            A_prev * (exp(-dt/tau1_prev) - exp(-dt/tau2_prev));
            tau1_prev;
            tau2_prev
        ];
        
        % 计算雅可比矩阵F
        F = zeros(3,3);
        F(1,1) = exp(-dt/tau1_prev) - exp(-dt/tau2_prev);
        F(1,2) = A_prev * (dt/tau1_prev^2) * exp(-dt/tau1_prev);
        F(1,3) = -A_prev * (dt/tau2_prev^2) * exp(-dt/tau2_prev);
        F(2,2) = 1;
        F(3,3) = 1;
        
        % 预测协方差
        P_pred = F * P * F' + Q;
        
        % ---------- 更新步骤 ----------
        % 观测模型：直接测量幅度
        H = [1, 0, 0];
        
        % 卡尔曼增益
        K = (P_pred * H') / (H * P_pred * H' + R);
        
        % 状态更新
        x_hat = x_pred + K * (raw_signal(k) - x_pred(1));
        P = (eye(3) - K * H) * P_pred;
        
        % 存储结果
        filtered_signal(k) = x_hat(1);
        x_estimates(:,k) = x_hat;
    end
end

