function filtered_signal = filter_xb(signal)
% 通用小波
    level = 4;          % 分解层数
    wavelet = 'db4';    % 小波基
    [c4, l4] = wavedec(signal, level, wavelet);
    sigma = median(abs(c4))/0.6745;
%     thr = sigma*sqrt(2*log(length(signal)));  % 通用阈值
    thr = sigma*sqrt(log(length(signal)));
    c_denoised4 = wthresh(c4, 's', thr);  % 使用小波软阈值
    filtered_signal = normalize(waverec(c_denoised4, l4, wavelet));


end
