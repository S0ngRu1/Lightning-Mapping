%函数：主辅窗口匹配
function [matched_windows ,tau_gcc12s]= find_match_windows(signal1, signal2, window_size, Fs)
    % 存储各个窗口的相关系数
    matched_windows = [];
    tau_gcc12s = [];
    window1 = signal1;
    for i = 1:256:2049
        % 提取当前窗口的数据
        window2 = signal2(i:i+window_size-1);
        [tau12,R12,lag12] = gccphat(window1,window2);
        max_R12 = maxvalue(R12);
        tau = lag12(maxindex(R12));
        tau_gcc12s = [tau_gcc12s; tau];
        matched_windows = [matched_windows; max_R12];
    end
end

