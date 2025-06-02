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

