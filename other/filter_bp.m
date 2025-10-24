
%% 设计巴特沃斯带通滤波器
function filtered_signal = filter_bp(signal, f1, f2, order)
    Fs = 400e6;          % 采样频率400MHz
    nyq = Fs / 2;        % 奈奎斯特频率=200MHz
    
    % 检查截止频率合理性（必须满足0 < f1 < f2 < nyq）
    if f1 <= 0 || f2 >= nyq || f1 >= f2
        error('截止频率错误：需满足 0 < f1 < f2 < %d MHz', nyq/1e6);
    end
    
    Wn = [f1, f2] / nyq; % 正确归一化：相对于奈奎斯特频率
    [b, a] = butter(order, Wn, 'bandpass'); % 明确指定带通类型（可选，但更清晰）
    filtered_signal = filtfilt(b, a, signal); % filtfilt零相位滤波，正确
end
