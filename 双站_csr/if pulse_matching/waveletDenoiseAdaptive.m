function denoised_signal = waveletDenoiseAdaptive(noisy_signal, level, wavelet)
% waveletDenoiseAdaptive - 改进的小波去噪（自适应阈值）
%
% 语法:
%   denoised_signal = waveletDenoiseAdaptive(noisy_signal, level, wavelet)
%
% 输入:
%   noisy_signal - 含噪信号向量
%   level        - 小波分解的层数
%   wavelet      - 小波基名称（如 'db4'）
%
% 输出:
%   denoised_signal - 去噪后重构的信号
%
% 示例:
%   level = 4;
%   wavelet = 'db4';
%   noisy_signal = chj3;  % 假设 chj3 已定义
%   denoised_signal = waveletDenoiseAdaptive(noisy_signal, level, wavelet);

    % 小波分解
    [c2, l2] = wavedec(noisy_signal, level, wavelet);
    c_denoised2 = c2;
    last = 0;
    
    % 对每一层细节系数进行处理
    for i = 1:level
        first = last + 1;
        last = first + l2(i) - 1;
        d = c2(first:last);
        
        % 自适应阈值计算
        sigma = median(abs(d)) / 0.6745;
        N = length(d);
        thr_univ = sigma * sqrt(2 * log(N));      % 通用阈值
        thr_mini = sigma * sqrt(2 * log(log(N)));   % minimax阈值
        
        % 混合阈值策略
        energy = sum(d.^2);
        if energy > median(d.^2) * N
            % 信号能量大，使用较小阈值
            thr_adapt = thr_mini;
        else
            % 信号能量小，使用较大阈值
            thr_adapt = thr_univ;
        end
        
        % 改进的阈值函数（平滑过渡）
        beta = 1.5;  % 平滑因子
        d_processed = sign(d) .* max(0, abs(d) - thr_adapt ./ (1 + exp(-beta*(abs(d)-thr_adapt))));
        
        c_denoised2(first:last) = d_processed;
    end
    
    % 小波重构得到去噪信号
    denoised_signal = waverec(c_denoised2, l2, wavelet);
end
