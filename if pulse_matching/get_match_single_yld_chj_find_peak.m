function [start_read_loc_chj, r_gccs] = get_match_single_yld_chj_find_peak(yld_signal_start_loc, skip_large_window)
    start_read_loc_chj = [];
    r_gccs = 0;
    window_lengths = [2e5, 2e4, 6000];

    if skip_large_window ~= 0
        window_lengths = window_lengths(2:end);
        current_chj_read_loc = skip_large_window;
    end

    % 对每个窗口长度进行匹配，逐步精细化
    for i = 1:length(window_lengths)
        current_window_length = window_lengths(i);
        % 读取yld信号
        yld_signal_length = current_window_length;
        yld_signal = read_signal('../20240822165932.6610CH1.dat', yld_signal_length, yld_signal_start_loc);
        filtered_yld_signal = filter_bp(yld_signal, 30e6, 80e6, 5);
        processed_yld_signal = real(windowsignal(detrend(filtered_yld_signal)));
        chj_length = current_window_length * 4;
        % 读取chj信号
        if i == 1 && skip_large_window == 0
            current_chj_read_loc = yld_signal_start_loc + 34371950 - current_window_length * 2;
            %             chj_length = current_window_length * 2;
        else
            current_chj_read_loc = current_chj_read_loc - current_window_length * 2;
        end
        if isempty(current_chj_read_loc)
            continue;
        end
        chj_signal = read_signal('../2024 822 85933.651462CH1.dat', chj_length, current_chj_read_loc);
        all_locs = [];
        all_R_gccs = [];
        all_t_gccs = [];
        filtered_chj_signal = filter_bp(chj_signal, 30e6, 80e6, 5);
        processed_chj_signal = real(windowsignal(detrend(filtered_chj_signal)));
        threshold = 30;
        % 寻找峰值
        [peaks, locs] = findpeaks(chj_signal, 'MinPeakHeight', threshold, 'MinPeakDistance', 256);
        all_locs = locs;
        % 遍历所有峰值
        num_peaks = numel(all_locs);
        if num_peaks == 0
            continue;
        end
        for pi = 1:num_peaks
            idx = all_locs(pi);
            % 确保峰值不超出信号范围
            if idx - (current_window_length / 2 - 1) <= 0 || idx + (current_window_length / 2) > length(processed_chj_signal)
                continue;
            end
            % 截取窗口信号
            signal1 = processed_chj_signal(idx - (current_window_length / 2 - 1):idx + (current_window_length / 2));
            % 去直流分量并应用窗函数
            ch1_new = real(windowsignal(detrend(signal1)));
            %互相关
            [r_gcc, lags_gcc] = xcorr(ch1_new, filtered_yld_signal, 'normalized');
            R_gcc = max(r_gcc);
            all_R_gccs = [all_R_gccs; R_gcc];
            t_gcc = cal_tau(r_gcc, lags_gcc');
            all_t_gccs = [all_t_gccs; t_gcc];
        end

        % 找到当前窗口中最大相关系数的位置
        [r_gccs, max_idx] = max(all_R_gccs);
        
        % 更新chj信号的起始位置：在当前读信号的位置上加上匹配得到的子窗口起始位置和时间偏移
        current_chj_start_loc = current_chj_read_loc + idx - current_window_length/2 + 1 + floor(all_t_gccs(max_idx));
        current_chj_read_loc = current_chj_start_loc;
        start_read_loc_chj = current_chj_start_loc;
    end
end