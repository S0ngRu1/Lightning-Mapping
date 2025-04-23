function [start_read_loc_chj, maxC] = get_match_single_yld_chj_find_peak(filtered_chj_signal1,filtered_yld_signal1,yld_signal_start_loc,r_loction)
    current_window_length = 1024;
%    绝对位置转相对位置
    yld_signal_start_loc = yld_signal_start_loc-r_loction;
    % 读取yld信号
    yld_signal_length = current_window_length;
    yld_signal = filtered_yld_signal1(yld_signal_start_loc+1-yld_signal_length/2: yld_signal_start_loc+yld_signal_length/2);
    processed_yld_signal = real(windowsignal(detrend(yld_signal)));
    loc_diff = 6000;
    chj_length = current_window_length + loc_diff;
    % 读取chj信号
    current_chj_read_loc = yld_signal_start_loc+1 - loc_diff/ 2;
    chj_signal = filtered_chj_signal1(current_chj_read_loc:chj_length+current_chj_read_loc);
    processed_chj_signal = real(windowsignal(detrend(chj_signal)));

    subplot(2,1,1);plot(processed_yld_signal);title('yld');xlabel('采样点数');ylabel('幅值');
    subplot(2,1,2);plot(processed_chj_signal);title('chj');xlabel('采样点数');ylabel('幅值');

    [r_gcc, lags_gcc] = xcorr(processed_yld_signal, processed_chj_signal,'none');
    t_gcc            = cal_tau(r_gcc, lags_gcc');
    r_gcc = r_gcc ./ (norm(processed_chj_signal) * norm(processed_yld_signal));
    [maxC, ~]   = max(r_gcc);
    start_read_loc_chj = current_chj_read_loc + t_gcc;
end