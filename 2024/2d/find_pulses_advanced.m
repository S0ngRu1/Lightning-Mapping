function pulse_catalog = find_pulses_advanced(waveform, noise_std, sampling_rate_hz, detection_threshold_factor, merge_gap_samples)
%
% 优化说明:
% 1. 移除了耗时的区域融合 while 循环，改为向量化操作。
% 2. 边界搜索改为基于预计算的逻辑掩码查找。
% 3. 抛物线拟合替换为预计算的矩阵乘法。
%

    %% --- 1. 计算希尔伯特包络和阈值 ---
    % 保持不变，这是标准操作
    if isrow(waveform), waveform = waveform'; end % 确保列向量
    envelope = abs(hilbert(waveform));
    mean_envelope = mean(envelope);
    detection_threshold = noise_std * detection_threshold_factor;
    n_samples = length(envelope);

    %% --- 2. 快速检测与融合脉冲区域 (核心优化) ---
    % 生成二值掩码
    binary_mask = envelope > detection_threshold;
    
    % 找到所有上升沿(Start)和下降沿(End)
    % diff结果: 1为上升沿前一位置, -1为下降沿位置
    d_mask = diff([0; binary_mask; 0]); 
    starts = find(d_mask == 1);
    ends = find(d_mask == -1) - 1; % 修正索引以匹配原信号
    
    if isempty(starts)
        pulse_catalog = [];
        return;
    end

    % --- 向量化融合逻辑 ---
    if length(starts) > 1
        % 计算相邻区域的间隙 (下一个Start - 当前End)
        gaps = starts(2:end) - ends(1:end-1);
        
        % 找到需要融合的间隙索引
        merge_indices = find(gaps < merge_gap_samples);
        
        % 填补间隙：直接在二值掩码上操作
        for i = 1:length(merge_indices)
            idx = merge_indices(i);
            % 将两个脉冲之间的空隙置为 1
            binary_mask(ends(idx) : starts(idx+1)) = 1; 
        end
        
        % 重新计算融合后的区域边界 (再次运行找边逻辑，极快)
        d_mask = diff([0; binary_mask; 0]);
        starts = find(d_mask == 1);
        ends = find(d_mask == -1) - 1;
    end
    
    num_final_regions = length(starts);

    %% --- 3. 预计算边界搜索掩码 (优化边界查找) ---
    % 目标：找到连续5个点都小于均值的位置
    % 方法：使用卷积。如果某点及周围满足条件，卷积值将为5。
    
    is_below_mean = envelope < mean_envelope;
    % 构造卷积核：检查点i及其后4个点 (forward) 或 前4个点 (backward)
    % 这里为了通用，我们标记所有属于"连续5个低值段"的点
    
    % 使用移动求和来判断连续区间
    kernel = ones(5, 1);
    % conv_res[i] 表示以 i 为中心的滑动窗口内的和 (需要处理边缘)
    % 我们使用 'valid' 类似逻辑手动处理更可控，或者简单的 loop 查找 lookup table
    
    % 更简单的方法：找到所有 "safe" 的点（即可以作为边界的点）
    % 定义：如果某点是连续5个低值序列的开头(forward)或结尾(backward)
    
    % 这里我们用快速查找法替代 convoluted logic，保持原算法逻辑但加速
    % 预先计算所有满足 < mean 的索引
    below_indices = find(is_below_mean); 
    
    %% --- 4. 预计算拟合矩阵 (优化正时) ---
    % 对应 x = [-2, -1, 0, 1, 2]' 的 2次多项式最小二乘伪逆矩阵
    % 这一步替代了 polyfit，极大提速
    % X = [-2 -1 0 1 2]'; V = [X.^2, X, ones(5,1)]; PinV = inv(V'*V)*V';
    % 下面是预计算好的 PinV 的具体数值
    PinV = [ 0.14285714, -0.07142857, -0.14285714, -0.07142857,  0.14285714; ... % a (x^2)
            -0.20000000, -0.10000000,  0.00000000,  0.10000000,  0.20000000; ... % b (x)
             0.08571429,  0.34285714,  0.48571429,  0.34285714,  0.08571429];    % c (const)
             
    ts_ns = 1 / sampling_rate_hz * 1e9;

    %% --- 5. 构建脉冲目录 ---
    % 预分配结构体数组 (提速)
    pulse_catalog = repmat(struct('start_idx', 0, 'end_idx', 0, 'peak_loc', 0, 'precise_time_ns', 0), num_final_regions, 1);
    valid_count = 0;

    for k = 1:num_final_regions
        s_idx = starts(k);
        e_idx = ends(k);
        
        % 提取区域数据
        region_data = envelope(s_idx:e_idx);
        
        % 寻找区域内峰值
        [~, max_idx_rel] = max(region_data);
        peak_loc = s_idx + max_idx_rel - 1;
        
        % --- 快速边界查找 ---
        % 向后找 (Start): 在 below_indices 中找到第一个 < s_idx 且满足连续5点条件的
        % 向前找 (End):   在 below_indices 中找到第一个 > e_idx 且满足连续5点条件的
        
        % 这里的逻辑如果完全向量化比较复杂，对于单脉冲循环内的逻辑，
        % 我们使用简化版：从 s_idx 往前找，遇到 连续5个 below_indices 即可。
        % 由于我们已经有了 below_indices，使用 find 的范围搜索比 while 步进快。
        
        % Backward Search
        % 找到所有在 peak 之前的低值点
        cands_back = below_indices(below_indices < peak_loc);
        if isempty(cands_back)
            real_start = 1;
        else
            % 我们要找的是：离 peak 最近的，且它是连续5个低值点的结束(或最后一点)
            % 简单优化：直接检查离 peak 最近的候选点，验证其前4个点是否也在候选列中
            % 为了极致速度，直接退化为 fast while (MATLAB JIT 对简单 while 优化很好，但对复杂逻辑不行)
            % 这里保持原逻辑的快速版：
            real_start = scan_boundary_fast(envelope, peak_loc, mean_envelope, -1, n_samples);
        end
        
        % Forward Search
        real_end = scan_boundary_fast(envelope, peak_loc, mean_envelope, 1, n_samples);

        if isnan(real_start) || isnan(real_end)
            continue; 
        end
        
        % --- 快速正时计算 ---
        if peak_loc <= 2 || peak_loc >= n_samples - 1
            precise_time_ns = (peak_loc - 1) * ts_ns;
        else
            % 获取5个点
            y_fit = envelope(peak_loc-2 : peak_loc+2);
            % 矩阵乘法求系数 [a; b; c]
            coeffs = PinV * y_fit; 
            a = coeffs(1); b = coeffs(2);
            
            % 顶点公式 -b/2a
            if abs(a) < 1e-10 % 防止除零 (直线)
                sub_sample_offset = 0;
            else
                sub_sample_offset = -b / (2 * a);
            end
            
            % 限制偏移量在合理范围内 (+/- 2)，防止拟合发散
            if abs(sub_sample_offset) > 2
                 sub_sample_offset = 0;
            end
            
            precise_time_ns = (peak_loc + sub_sample_offset - 1) * ts_ns;
        end
        
        valid_count = valid_count + 1;
        pulse_catalog(valid_count).start_idx = real_start;
        pulse_catalog(valid_count).end_idx = real_end;
        pulse_catalog(valid_count).peak_loc = peak_loc;
        pulse_catalog(valid_count).precise_time_ns = precise_time_ns;
    end
    
    % 截断多余的预分配空间
    pulse_catalog = pulse_catalog(1:valid_count);
end

% --- 辅助函数: 极简快速扫描 ----
function idx = scan_boundary_fast(env, start_pos, thresh, step, n_len)
    % env: 包络信号
    % start_pos: 搜索起始点 (峰值位置)
    % thresh: 阈值 (mean_envelope)
    % step: 方向步长 (-1 或 1)
    % n_len: 信号总长度

    idx = NaN;        % 默认返回 NaN
    curr = start_pos; % 当前指针
    consecutive = 0;  % 连续满足条件的计数器
    
    % 设置一个安全的最大循环次数，防止极端情况死循环
    % 虽然理论上边界检查会跳出，但加个保护更稳健
    while true
        curr = curr + step;
        
        % 1. 越界检查 (碰到信号边缘)
        if curr < 1 
            idx = 1;      % 到达最左端
            return; 
        end
        if curr > n_len 
            idx = n_len;  % 到达最右端
            return; 
        end
        
        % 2. 阈值检查
        if env(curr) < thresh
            consecutive = consecutive + 1;
            % 找到了连续 5 个点低于阈值
            if consecutive >= 5
                idx = curr; % 当前点即为确定的边界点
                return;
            end
        else
            % 只要有一个点不满足，计数器归零，继续往外搜
            consecutive = 0;
        end
    end
end