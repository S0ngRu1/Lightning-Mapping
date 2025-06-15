function windowed_signal = windowsignal(signal)
%     r_length = length(signal);
%    % 使用汉明窗
%    window = hamming(r_length);
%    % 对滤波后的信号应用窗函数
%    windowed_signal = signal .* window; % 信号与窗函数相乘

    X = fft(signal);      %变换到频域加窗
    r_length = length(X);
    window = hamming(r_length);
%     得到的是频域信号
    X_windowed = X .* window;

% %     % 进行逆傅里叶变换得到时域信号
      windowed_signal = ifft(X_windowed);

end

