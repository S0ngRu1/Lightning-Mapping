function [start_read_loc_chj, r_gccs] = get_match_single_yld_chj_find_peak(filtered_chj_signal1,filtered_yld_signal1,yld_signal_start_loc, skip_large_window)
    start_read_loc_chj = [];
    r_gccs = 0;
    window_lengths = [2e5, 2e4, 6000];%200us,20us,6us

    if skip_large_window ~= 0
        window_lengths = window_lengths(2:end);
        current_chj_read_loc = skip_large_window;
    end
    yld_signal_start_loc = yld_signal_start_loc-3e8-1.5e4 +1;
    % 对每个窗口长度进行匹配，逐步精细化
    for i = 1:length(window_lengths)
        current_window_length = window_lengths(i);
        % 读取yld信号
        yld_signal_length = current_window_length;
        
        processed_yld_signal = filtered_yld_signal1(yld_signal_start_loc: yld_signal_start_loc+yld_signal_length-1);
        chj_length = current_window_length * 4;
        % 读取chj信号
        if i == 1 && skip_large_window == 0
            current_chj_read_loc = yld_signal_start_loc - current_window_length * 2;
            %             chj_length = current_window_length * 2;
        else
            current_chj_read_loc = current_chj_read_loc - current_window_length * 2;
        end
        if isempty(current_chj_read_loc)
            continue;
        end
        chj_signal = filtered_chj_signal1(current_chj_read_loc:current_chj_read_loc + chj_length-1);
        all_locs = [];
        all_R_gccs = [];
        all_t_gccs = [];
        % 寻峰匹配计算所有可能的三维源
        noise = read_signal('../2024 822 85933.651462CH1.dat',65400,1e8);
        filtered_chj_noise = rfi_filter(noise,65400);
        threshold = 5*std(filtered_chj_noise);
      
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
            if idx - current_window_length / 2  <= 0 || idx + (current_window_length / 2) -1 > length(chj_signal)
                continue;
            end
            % 截取窗口信号
            ch1_new = chj_signal(idx - (current_window_length / 2 ):idx + (current_window_length / 2)- 1);
            %互相关
            [r_gcc, lags_gcc] = xcorr(ch1_new, processed_yld_signal, 'normalized');
            R_gcc = max(r_gcc);
            all_R_gccs = [all_R_gccs; R_gcc];
            t_gcc = cal_tau(r_gcc, lags_gcc');
            all_t_gccs = [all_t_gccs; t_gcc];
        end

        % 找到当前窗口中最大相关系数的位置
        [r_gccs, max_idx] = max(all_R_gccs);
        
        % 更新chj信号的起始位置：在当前读信号的位置上加上匹配得到的子窗口起始位置和时间偏移
        current_chj_start_loc = current_chj_read_loc + idx - current_window_length/2 + floor(all_t_gccs(max_idx));
        current_chj_read_loc = current_chj_start_loc;
        start_read_loc_chj = current_chj_start_loc;
    end
end