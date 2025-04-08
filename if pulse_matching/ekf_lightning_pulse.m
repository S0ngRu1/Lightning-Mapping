function [filtered_signal, x_estimates] = ekf_lightning_pulse(raw_signal, fs, x0, Q, R, P0)
    % EKF_LIGHTNING_PULSE 扩展卡尔曼滤波处理闪电脉冲信号
    %
    % 输入参数：
    %   raw_signal : 原始信号向量 (Nx1)
    %   fs         : 采样频率 (Hz)
    %   x0         : 初始状态 [A0; τ1] (2x1)
    %   Q          : 过程噪声协方差矩阵 (2x2)
    %   R          : 观测噪声协方差 (标量)
    %   P0         : 初始估计协方差矩阵 (2x2)
    %
    % 输出参数：
    %   filtered_signal : 滤波后信号 (Nx1)
    %   x_estimates     : 状态估计轨迹 (2xN)

    dt = 1/fs;              % 采样间隔
    n_samples = length(raw_signal);
    filtered_signal = zeros(n_samples, 1);
    x_estimates = zeros(2, n_samples);
    
    % 初始化
    x_hat = x0(:);          % 确保列向量
    P = P0;
    
    for k = 1:n_samples
        % --- 预测步骤 ---
        A_prev = x_hat(1);
        tau1_prev = x_hat(2);
        
        % 非线性状态预测（双指数衰减模型）
        x_pred = [A_prev * exp(-dt/tau1_prev); 
                  tau1_prev];
        
        % 计算雅可比矩阵F
        F = [exp(-dt/tau1_prev),    A_prev*(dt/tau1_prev^2)*exp(-dt/tau1_prev);
             0,                      1];
        
        % 预测协方差
        P_pred = F * P * F' + Q;
        
        % --- 更新步骤 ---
        % 观测预测（假设直接观测幅度）
        z_pred = x_pred(1);
        
        % 观测雅可比矩阵H
        H = [1, 0];
        
        % 卡尔曼增益
        K = P_pred * H' / (H * P_pred * H' + R);
        
        % 状态更新
        x_hat = x_pred + K * (raw_signal(k) - z_pred);
        P = (eye(2) - K * H) * P_pred;
        
        % 存储结果
        filtered_signal(k) = x_hat(1);
        x_estimates(:, k) = x_hat;
    end
end