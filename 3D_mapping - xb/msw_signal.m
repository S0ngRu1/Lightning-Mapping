function mswed_signal = msw_signal(signal , peak_x ,msw_length,signal_length)
      % 找到峰值的 x 值在信号中的索引
    left_idx = max(peak_x - msw_length+1, 1);  % 确定左边界的索引
    right_idx = min(peak_x + msw_length, signal_length);  % 确定右边界的索引
    mswed_signal = signal(left_idx:right_idx);  % 提取以中心 x 值为中心的左右40个采样点

end