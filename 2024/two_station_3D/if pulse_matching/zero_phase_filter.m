function x_filtered = zero_phase_filter(x, fs, corner_frequency)
    % 零相位滤波器
    % x: 输入信号
    % fs: 采样频率 (Hz)
    % corner_frequency: 角频率 (Hz)
    % x_filtered: 滤波后的信号
    
    Ts = 1/fs; % 采样周期 (秒)
    
    % 计算滤波器系数
    alpha = tan(pi * corner_frequency * Ts);
    b0 = 1 / (1 + alpha);
    b1 = 2 * b0;
    b2 = b0;
    a1 = -2 * b0;
    a2 = b0;
    
    % 设计单级滤波器
    b = [b0, b1, b2];
    a = [1, a1, a2];
    
    % 应用零相位滤波
    x_filtered = filtfilt(b, a, x);
end
