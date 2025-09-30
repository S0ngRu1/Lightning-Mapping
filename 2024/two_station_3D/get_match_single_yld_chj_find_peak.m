function [start_read_loc_chj, maxC] = get_match_single_yld_chj_find_peak(filtered_chj_signal1,filtered_yld_signal1,yld_signal_start_loc,r_loction)
    current_window_length = 10240;
%    绝对位置转相对位置
    yld_signal_start_loc = yld_signal_start_loc-r_loction;
    % 读取yld信号
    yld_signal_length = current_window_length;
    yld_signal = filtered_yld_signal1(yld_signal_start_loc+1-yld_signal_length/2: yld_signal_start_loc+yld_signal_length/2);%yld_signal_start_loc为窗口起始位置，-yld_signal_length/2问了将辐射源脉冲置于窗口中间
    processed_yld_signal = real(windowsignal(detrend(yld_signal)));
    loc_diff = 40000;
    chj_length = current_window_length + loc_diff;
    % 读取chj信号
    current_chj_read_loc = yld_signal_start_loc+1 - loc_diff/ 2;%在与引雷点对应位置取10240，并前后多取20000点作为遍历范围，在其中找匹配窗口
    chj_signal = filtered_chj_signal1(current_chj_read_loc:current_chj_read_loc+chj_length);
    processed_chj_signal = real(windowsignal(detrend(chj_signal)));
    [r_gcc, lags_gcc] = xcorr(processed_yld_signal, processed_chj_signal,'none');
    t_gcc            = cal_tau(r_gcc, lags_gcc');
    r_gcc = r_gcc ./ (norm(processed_chj_signal) * norm(processed_yld_signal));
    [maxC, ~]   = max(r_gcc);
    start_read_loc_chj = current_chj_read_loc + t_gcc;
end