function mswed_signal = msw_signal(signal , peak_x ,length)
    peak_x_idx = find(signal(:,1) == peak_x);  % 找到峰值的 x 值在信号中的索引
    left_idx = max(peak_x_idx - length+1, 1);  % 确定左边界的索引
    right_idx = min(peak_x_idx + length, 10240);  % 确定右边界的索引
    mswed_signal = signal(left_idx:right_idx,2);  % 提取以中心 x 值为中心的左右40个采样点

end