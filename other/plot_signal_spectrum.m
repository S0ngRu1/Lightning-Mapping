function plot_signal_spectrum(signal)
% 绘制时域信号的频谱图（修正版）
% 输入：
%   signal - 输入时域信号
    fs = 400e6;
    n = length(signal);
    % 计算FFT（补零至下一个2的幂次，提高频谱分辨率，可选）
    fft_signal = fft(signal);
    
    % 1. 计算频率轴（0 ~ fs/2，正频率范围）
    f = (0:n/2) * (fs / n);  % 频率间隔为 fs/n，共 n/2+1 个点（含0和fs/2）
    
    % 2. 计算振幅（归一化）
    amp = abs(fft_signal(1:n/2+1));  % 取正频率部分
    amp(2:end-1) = 2 * amp(2:end-1) / n;  % 非直流/奈奎斯特频率点：乘以2/n
    amp(1) = amp(1) / n;  % 直流分量：仅除以n（无负频率对应）
    if mod(n, 2) == 0  % 偶数点时，奈奎斯特频率点无需加倍
        amp(end) = amp(end) / n;
    end
    
    % 3. 绘图
    figure;
    plot(f/1e6, amp);  % 若fs单位为Hz，转换为MHz显示（除以1e6）
    xlabel('Frequency (MHz)');
    ylabel('Amplitude');
    grid on;
    title('Signal Spectrum');
end