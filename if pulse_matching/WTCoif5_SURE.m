function denoised = WTCoif5_SURE(signal, level)
    % 使用 Coif5 小波和 SURE 阈值去噪
    % signal: 输入信号（1D）
    % level: 分解层数（建议为 log2(length(signal)) 的一部分）

    % 设置小波基
    wavelet = 'coif5';
    
    % 小波分解
    [C, L] = wavedec(signal, level, wavelet);
    
    % 获取细节系数并估计噪声标准差
    detail_coeffs = detcoef(C, L, level);  % 最底层细节
    sigma = median(abs(detail_coeffs)) / 0.6745;  % 鲁棒估计
    
    % 初始化去噪系数
    C_denoised = C;
    
    % 对每一层细节系数进行 SURE 去噪
    for j = 1:level
        Dj = detcoef(C, L, j);
        idx = sum(L(1:end-j));  % 找到 Dj 在 C 中的起始位置
        Dj_len = length(Dj);
        pos = idx + (1:Dj_len);
        
        % 计算 SURE 阈值
        thr = sure_threshold(Dj, sigma);
        
        % 应用软阈值
        C_denoised(pos) = wthresh(Dj, 's', thr);
    end
    
    % 信号重构
    denoised = waverec(C_denoised, L, wavelet);
end

function thr = sure_threshold(coeff, sigma)
    % SURE 阈值估计 (软阈值)
    n = length(coeff);
    coeff2 = sort(abs(coeff).^2);
    risks = zeros(n,1);
    for i = 1:n
        t = sqrt(coeff2(i));
        risks(i) = n - 2*i + sum(min(coeff2, coeff2(i))) + i*coeff2(i);
    end
    [~, idx] = min(risks);
    thr = sqrt(coeff2(idx)) * sigma;
end
