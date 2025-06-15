function [gcc, lags] = gcc_phat(x, y, fs)
    % 广义互相关 - PHAT 版本
    % x, y: 输入信号（可以是不同长度）
    % fs: 采样率（用于计算真实时间延迟，可选）
    
    N = length(x) + length(y);  % 确定变换长度（零填充后避免圆周卷积）
    
    % FFT
    X = fft(x, N);
    Y = fft(y, N);
    
    % 互功率谱
    Gxy = X .* conj(Y);
    
    % PHAT 加权
    Gxy_phat = Gxy ./ (abs(Gxy) + eps);  % 防止除以零
    
    % IFFT 得到互相关函数
    gcc = real(ifft(Gxy_phat));
    
    % 使输出与延迟对齐
    gcc = fftshift(gcc);
    lags = (-floor(N/2):ceil(N/2)-1);  % 延迟索引
    
    % 若提供采样率，则将 lags 转换为秒
    if nargin == 3
        lags = lags / fs;
    end
end
