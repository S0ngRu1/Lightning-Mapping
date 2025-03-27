function [start_read_loc_chj, r_gccs] = get_match_single_yld_chj_siddle(yld_signal_start_loc, skip_large_window)
    % 定义逐步缩小的窗口长度，依次进行粗匹配到细匹配
    % 小波去噪参数
%     level = 7;          % 分解层数
%     wavelet = 'db4';    % 小波基 
    window_lengths = [2e5, 2e4, 6000];
    % 如果skip_large_window参数不为0，则跳过第一个窗口长度（2e7）
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
        filtered_yld_signal = filter_bp(yld_signal, 20e6, 80e6, 5);
        processed_yld_signal = real(windowsignal(detrend(filtered_yld_signal)));
        chj_length = current_window_length * 4;
        % 读取chj信号
        if i == 1 && skip_large_window == 0
            current_chj_read_loc = yld_signal_start_loc + 34371950 - current_window_length * 2;
%             chj_length = current_window_length * 2;
        else
            current_chj_read_loc = current_chj_read_loc - current_window_length * 2;
        end
        chj_signal = read_signal('../2024 822 85933.651462CH1.dat', chj_length, current_chj_read_loc);
        % 设置滑动窗口参数（步长设为yld信号长度的1/4）
        subsignal_step = yld_signal_length / 10;
        subsignal_starts = 1:subsignal_step:chj_length;
        all_R_gccs = [];
        all_t_gccs = [];
        % 对当前窗口内的每个子段进行匹配
        for subi = 1:numel(subsignal_starts)
            if subsignal_starts(subi) + yld_signal_length - 1 > chj_length
                continue
            end
            subsignal_chj = chj_signal(subsignal_starts(subi) : subsignal_starts(subi) + yld_signal_length - 1);
            filtered_chj_signal = filter_bp(subsignal_chj, 20e6, 80e6, 5);
%             filtered_chj_signal = waveletDenoiseAdaptive(subsignal_chj, level, wavelet);
            processed_chj_signal = real(windowsignal(detrend(filtered_chj_signal)));
            
            [r_gcc, lags_gcc] = xcorr(processed_chj_signal, processed_yld_signal, 'normalized');
            R_gcc = max(r_gcc);
            all_R_gccs = [all_R_gccs; R_gcc];
            t_gcc = cal_tau(r_gcc, lags_gcc');
            all_t_gccs = [all_t_gccs; t_gcc];
        end

        % 找到当前窗口中最大相关系数的位置
        [r_gccs, max_idx] = max(all_R_gccs);
        [r_gccsk, max_idxk] = maxk(all_R_gccs,5);
        % 更新chj信号的起始位置：在当前读信号的位置上加上匹配得到的子窗口起始位置和时间偏移
        current_chj_start_loc = current_chj_read_loc + subsignal_starts(max_idx) + floor(all_t_gccs(max_idx));
        current_chj_read_loc = current_chj_start_loc;
    end
    % 最终匹配位置和最后一次的相关系数作为输出
    start_read_loc_chj = current_chj_start_loc;
end
