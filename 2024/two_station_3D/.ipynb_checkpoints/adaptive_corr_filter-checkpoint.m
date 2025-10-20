function logicalIndex = adaptive_corr_filter(resultTable, winLenCol, corrCol, baseQuantile)
    % 离散长度映射：针对每种长度分别筛选高相关系数
    % resultTable: 输入表
    % winLenCol: 窗口长度列名（如 'Length'）
    % corrCol: 相关系数列名（如 'Rcorr'）
    % baseQuantile: 分位数阈值（如 0.8）

    lengths = resultTable.(winLenCol);
    rcorrs = resultTable.(corrCol);
    logicalIndex = false(height(resultTable), 1);

    unique_lengths = unique(lengths);

    for i = 1:length(unique_lengths)
        len = unique_lengths(i);
        bin_mask = lengths == len;
        bin_rcorr = rcorrs(bin_mask);

        if length(bin_rcorr) >= 3  % 至少3个才能用分位数比较稳定
            threshold = quantile(bin_rcorr, baseQuantile);
            logicalIndex(bin_mask) = bin_rcorr >= threshold;
        elseif length(bin_rcorr) == 2
            logicalIndex(bin_mask) = bin_rcorr >= min(bin_rcorr);  % 取大的
        elseif length(bin_rcorr) == 1
            logicalIndex(bin_mask) = true;  % 只有一个，保留
        end
    end
end
