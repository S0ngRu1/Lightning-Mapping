function [start_read_loc_chj, r_gccs] = get_match_single_yld_chj_find_peak(filtered_chj_signal1,filtered_yld_signal1,yld_signal_start_loc, skip_large_window)
    % 定义逐步缩小的窗口长度，依次进行粗匹配到细匹配
    start_read_loc_chj = [];
    r_gccs = 0;    
    window_lengths = [ 2e5, 2e4,  6000,1024];
    % 如果skip_large_window参数不为0，则跳过第一个窗口长度（2e7）
    if skip_large_window ~= 0
        window_lengths = window_lengths(2:end);
        current_chj_read_loc = skip_large_window;
    end
%    绝对位置转相对位置
    yld_signal_start_loc = yld_signal_start_loc-4.4e8;
    % 对每个窗口长度进行匹配，逐步精细化
    for i = 1:length(window_lengths)
        current_window_length = window_lengths(i);
        % 读取yld信号
        yld_signal_length = current_window_length;
        processed_yld_signal = filtered_yld_signal1(yld_signal_start_loc+1: yld_signal_start_loc+yld_signal_length);
        processed_yld_signal = real(windowsignal(detrend(processed_yld_signal)));
        loc_diff = min(2*1.5e5,2*1.5e5/(2e5/current_window_length));
        chj_length = current_window_length + loc_diff;
        % 读取chj信号
        if i == 1 && skip_large_window == 0
            current_chj_read_loc = yld_signal_start_loc - loc_diff/ 2;
            %             chj_length = current_window_length * 2;
        else
            current_chj_read_loc = current_chj_read_loc - loc_diff/ 2;
        end
        if isempty(current_chj_read_loc)
            continue;
        end
        if current_chj_read_loc<0 || chj_length+current_chj_read_loc > length(filtered_chj_signal1)
            continue
        end
        chj_signal = filtered_chj_signal1(current_chj_read_loc+1:chj_length+current_chj_read_loc);
        [r_gcc, lags_gcc] = gcc_phat(normalize(chj_signal), normalize(processed_yld_signal));
        t_gcc = cal_tau(r_gcc, lags_gcc');
        % 更新chj信号的起始位置
        current_chj_start_loc = current_chj_read_loc - floor(t_gcc);
        current_chj_read_loc = current_chj_start_loc;
        start_read_loc_chj = current_chj_start_loc;
    end
end