function [filtered_signal, dominant_freqs] = extract_band_emd(signal, fs, freq_lower, freq_upper, plot_flag)
% extract_band_emd: 从EMD中提取主导频率在指定频段内的IMFs并重构信号
%
% 输入：
%   signal       - 一维输入信号（列向量）
%   fs           - 采样频率 (Hz)，例如 200e6
%   freq_lower   - 主导频率下限 (Hz)，例如 10e6
%   freq_upper   - 主导频率上限 (Hz)，例如 100e6
%   plot_flag    - 是否绘图（true/false）
%
% 输出：
%   filtered_signal - 重构的带通滤波信号
%   dominant_freqs  - 每个IMF的主导频率（Hz）

    signal = signal(:);
    [imfs, ~] = emd(signal, 'Interpolation', 'spline');

    num_imfs = size(imfs, 2);
    dominant_freqs = zeros(1, num_imfs);
    selected_imfs = [];

    for i = 1:num_imfs
        imf = imfs(:, i);
        N = length(imf);
        Y = fft(imf);
        f = (0:N-1)*(fs/N);

        half_N = floor(N/2);
        magnitude = abs(Y(1:half_N));
        freq_axis = f(1:half_N);

        % 主导频率 = 最大幅值对应频率
        [~, max_idx] = max(magnitude);
        dominant_freq = freq_axis(max_idx);
        dominant_freqs(i) = dominant_freq;

        % 判断是否落在频段内
        if dominant_freq >= freq_lower && dominant_freq <= freq_upper
            selected_imfs(:, end+1) = imf;
        end
    end

    % 重构带通信号
    if isempty(selected_imfs)
        filtered_signal = zeros(size(signal));
        warning('没有主导频率落在 %.2f MHz ~ %.2f MHz 范围内的IMF。', ...
            freq_lower/1e6, freq_upper/1e6);
    else
        filtered_signal = sum(selected_imfs, 2);
    end

    % 可视化
    if plot_flag
        figure('Name', 'IMF Frequency Spectrum', 'NumberTitle', 'off');
        for i = 1:num_imfs
            subplot(ceil(num_imfs/2), 2, i);
            N = length(imfs(:, i));
            f = (0:N-1)*(fs/N);
            half_N = floor(N/2);
            magnitude = abs(fft(imfs(:, i)));
            plot(f(1:half_N)/1e6, magnitude(1:half_N));
            hold on;
            xline(dominant_freqs(i)/1e6, 'r--', sprintf('%.1f MHz', dominant_freqs(i)/1e6));
            xline(freq_lower/1e6, 'k:', 'Low');
            xline(freq_upper/1e6, 'k:', 'High');
            title(sprintf('IMF %d', i));
            xlabel('Frequency (MHz)');
            ylabel('Magnitude');
            grid on;
        end
        sgtitle('IMF 主导频率筛选可视化');
    end
end
