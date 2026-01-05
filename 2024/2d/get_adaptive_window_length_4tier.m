function window_len = get_adaptive_window_length_4tier(pulse_count,single_length)
    thresholds = [60, 120, 250] * single_length/1e4;
    T1 = thresholds(1);
    T2 = thresholds(2);
    T3 = thresholds(3);
    
    if pulse_count <= T1
        window_len = 4096;
    elseif pulse_count > T1 && pulse_count <= T2
        window_len = 2048;
    elseif pulse_count > T2 && pulse_count <= T3
        window_len = 1024;
    else % pulse_count > T3
        window_len = 512;
    end
end
