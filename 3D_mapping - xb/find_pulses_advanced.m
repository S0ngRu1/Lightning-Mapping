function pulse_catalog = find_pulses_advanced(waveform, noise_std, sampling_rate_hz, detection_threshold_factor,merge_gap_samples)
%
% 输入:
%   waveform        - 需要进行脉冲发现的信号波形
%   noise_std        - 噪声水平
%   sampling_rate_hz - 信号的采样率 (Hz)
%   detection_threshold_factor - 包络检测阈值因
%
% 输出:
%   pulse_catalog   - 结构体数组, 每个元素代表一个脉冲，包含:
%                     .start_idx      - 脉冲的起始采样点索引
%                     .end_idx        - 脉冲的结束采样点索引
%                     .peak_loc       - 脉冲内部最高峰的索引
%                     .precise_time_ns - 抛物线拟合得到的亚采样点精度时刻 (ns)
    
     %% --- 1. 计算希尔伯特包络和阈值 ---
    envelope = abs(hilbert(waveform));
    mean_envelope = mean(envelope);
    detection_threshold = noise_std * detection_threshold_factor;
    
    %% --- 2. 初步检测脉冲区域 ---
    pulse_regions = bwlabel(envelope > detection_threshold);
    num_initial_regions = max(pulse_regions);
    
    if num_initial_regions < 2
        regions_to_process = pulse_regions;
    else
        % ************* 融合邻近的脉冲区域 *************
        fprintf('初步检测到 %d 个区域，开始进行邻近区域融合...\n', num_initial_regions);
        
        merged_regions = pulse_regions; 
        
        while true % 持续循环，直到某一轮没有任何融合发生
            
            was_merged_in_this_pass = false; % 本轮融合的标志
            
            unique_labels = unique(merged_regions);
            unique_labels = unique_labels(unique_labels > 0); % 去掉背景0
            
            if numel(unique_labels) < 2
                break; % 如果只剩一个或没有区域，则结束
            end

            % 获取每个区域的边界
            stats = regionprops(merged_regions', 'BoundingBox');
            region_boundaries = zeros(numel(stats), 2);
            for k = 1:numel(stats)
                bbox = stats(k).BoundingBox;
                % --- 修正1：边界计算的“差一”错误 ---
                region_boundaries(k, 1) = floor(bbox(1));
                region_boundaries(k, 2) = floor(bbox(1)) + floor(bbox(3)) - 1; % 正确的终点
            end
            
            % 检查相邻区域的间隙
            for k = 1:numel(unique_labels) - 1
                current_label = unique_labels(k);
                next_label = unique_labels(k+1);
                
                % 找到属于这两个标签的所有点的索引
                current_indices = find(merged_regions == current_label);
                next_indices = find(merged_regions == next_label);
                
                current_region_end = max(current_indices);
                next_region_start = min(next_indices);
                
                gap = next_region_start - current_region_end;
                
                if gap < merge_gap_samples
                    % 将下一个区域合并到当前区域
                    merged_regions(merged_regions == next_label) = current_label;
                    was_merged_in_this_pass = true;
                    break; % 发生了一次融合，立即跳出内层循环，重新开始外层while循环
                end
            end
            
            if ~was_merged_in_this_pass
                % 如果完整的一轮扫描都没有发生任何融合，则说明融合已完成
                break;
            end
        end % while true
        
        % 对融合后的区域重新编号，使其连续
        [~, ~, regions_to_process] = unique(merged_regions);
        regions_to_process = reshape(regions_to_process, size(pulse_regions));
        regions_to_process(pulse_regions==0) = 0; % 恢复背景
    end
    
    num_final_regions = max(regions_to_process);
    
    %% --- 3. 对最终的脉冲区域进行精确界定和正时 (不变) ---
    pulse_catalog = [];
    
    for k = 1:num_final_regions
        region_indices = find(regions_to_process == k);
        if isempty(region_indices), continue; end
        [~, max_idx_in_region] = max(envelope(region_indices));
        peak_loc = region_indices(max_idx_in_region);
        start_idx = find_pulse_boundary(envelope, region_indices(1), mean_envelope, 'backward');
        end_idx = find_pulse_boundary(envelope, region_indices(end), mean_envelope, 'forward');
        if isnan(start_idx) || isnan(end_idx), continue; end
        precise_time_ns = get_precise_timing(envelope, peak_loc, sampling_rate_hz);
        pulse.start_idx = start_idx;
        pulse.end_idx = end_idx;
        pulse.peak_loc = peak_loc;
        pulse.precise_time_ns = precise_time_ns;
        if isempty(pulse_catalog), pulse_catalog = pulse; else, pulse_catalog(end+1) = pulse; end
    end
end

% --- 辅助函数1: 寻找脉冲边界 ---
function boundary_idx = find_pulse_boundary(envelope, peak_loc, mean_envelope, direction)
    search_loc = peak_loc;
    boundary_idx = NaN;
    max_iter = length(envelope); % 防止无限循环
    iter_count = 0;

    while iter_count < max_iter
        if strcmpi(direction, 'backward')
            % 检查峰值之前的5个点
            check_indices = (search_loc - 5) : (search_loc - 1);
            if check_indices(1) < 1
                boundary_idx = 1; % 到达信号开头
                break;
            end
            
            if all(envelope(check_indices) < mean_envelope)
                boundary_idx = check_indices(1);
                break;
            else
                search_loc = search_loc - 1;
            end
        else % forward
            % 检查峰值之后的5个点
            check_indices = (search_loc + 1) : (search_loc + 5);
            if check_indices(end) > length(envelope)
                boundary_idx = length(envelope); % 到达信号结尾
                break;
            end
            
            if all(envelope(check_indices) < mean_envelope)
                boundary_idx = check_indices(end);
                break;
            else
                search_loc = search_loc + 1;
            end
        end
        iter_count = iter_count + 1;
    end
end

% --- 辅助函数2: 抛物线拟合正时 ---
function precise_time_ns = get_precise_timing(envelope, peak_loc, sampling_rate_hz)
    ts_ns = 1 / sampling_rate_hz * 1e9;

    % 检查边界，确保可以取到5个点
    if peak_loc <= 2 || peak_loc >= length(envelope) - 1
        % 如果点太靠边，无法取5个点，则直接返回峰值点的时刻
        precise_time_ns = (peak_loc - 1) * ts_ns;
        return;
    end
    
    % 取以最高峰为中心的5个包络点
    y_fit = envelope(peak_loc-2 : peak_loc+2);
    x_fit = (-2:2)';
    
    % 进行二次多项式(抛物线)拟合
    p_coeffs = polyfit(x_fit, y_fit, 2);
    
    % 抛物线顶点公式 x = -b / (2a)
    % p_coeffs = [a, b, c]
    sub_sample_offset = -p_coeffs(2) / (2 * p_coeffs(1));
    
    % 计算亚采样点精度的位置
    precise_loc = peak_loc + sub_sample_offset;
    
    % 转换为纳秒
    precise_time_ns = (precise_loc - 1) * ts_ns;
end