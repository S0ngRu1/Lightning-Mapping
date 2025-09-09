function delay = cal_delay(R_xy)
    r_xy = ifft(R_xy);

% 找到主峰值位置
[~, max_idx] = max(r_xy);

% 计算估计的时间延迟
delay = max_idx / 200e3;

end