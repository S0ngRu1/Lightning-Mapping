function filtered_signal = rfi_filter(ori_signal,sub_signal_length)
    filtered_signal = [];
    subsignal_starts = 1:sub_signal_length/2:length(ori_signal);
    for i = 1:length(subsignal_starts)
        if subsignal_starts(i) + sub_signal_length - 1 > length(ori_signal)
            continue
        end
        subsignal = ori_signal(subsignal_starts(i):subsignal_starts(i)+sub_signal_length-1);
        windowed_signal = window_plus(sub_signal_length,subsignal);
        sub_filtered_signal = datafilter(windowed_signal);
        filtered_signal = [filtered_signal; sub_filtered_signal(sub_signal_length*0.25+1:sub_signal_length*0.75)];
    end
end