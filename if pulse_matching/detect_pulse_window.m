
function [pulse_flag] = detect_pulse_window(signal)
    N = round(200e6 * 1e-6); 
    P = abs(hilbert(signal)).^2; % 网页8的功率计算
    P_smooth = movmean(P, N); % 滑动平均
    
    % 双窗口能量比
    for k = 2*N:length(P_smooth)
        P_A = mean(P_smooth(k-2*N+1:k-N));
        P_B = mean(P_smooth(k-N:k));
        ratio_up(k) = P_B / (P_A + eps);
        ratio_down(k) = P_A / (P_B + eps);
    end
    
    % 阈值判决
    rise_edges = find(ratio_up > 1.5, 1);
    fall_edges = find(ratio_down > 1.3, 1);
    pulse_flag = ~isempty(rise_edges) && ~isempty(fall_edges);
end