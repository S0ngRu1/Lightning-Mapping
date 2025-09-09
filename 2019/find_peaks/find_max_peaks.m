%函数：寻找信号的最大峰的位置
function max_peak_loc = find_max_peaks(x)
    % 找到信号的峰值
    [peaks, peak_locs] = findpeaks(x);
    [~, max_idx] = max(peaks);
    max_peak_loc = peak_locs(max_idx);
end