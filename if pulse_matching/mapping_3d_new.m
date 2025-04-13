%% Phase 0: 参数初始化 ============================================
fprintf('Phase 0: 参数初始化\n');
clear; clc; close all;
Fs = 200e6; 
yld_data_path = '..\\';
chj_data_path = '..\\'; 
yld_base_filename = '20240822165932.6610'; 
chj_base_filename = '2024 822 85933.651462';

% --- 站点信息 ---
yld_sit = [0, 0, 0];            
chj_sit = [1991, -7841.2, 0];
dist_m = norm(chj_sit - yld_sit); 
c_mpns = 0.299792458;             
delta_t_theo_ns = dist_m / c_mpns; % 理论传播时间
delta_t_theo_samples = round(delta_t_theo_ns * Fs/5); % 理论传播时间（样本数）

% --- 粗略对齐参数 ---
coarse_align_yld_ch_idx = 1; % YLD 用于粗略对齐的天线
coarse_align_chj_ch_idx = 1; % CHJ 用于粗略对齐的天线
% ** 通过绘图得到近似的位置 **
anchor_loc_yld_approx = 4.8e8; % YLD 脉冲事件的近似样本索引
anchor_loc_chj_approx = 4.8e8 + 3.4e7; % CHJ 脉冲事件的近似样本索引
W_coarse_duration = 200e-6; % YLD 脉冲窗口的持续时间
W_search_duration = 200e-3; % CHJ 搜索窗口的持续时间（设置大一些）
W_coarse_half_samples = round(W_coarse_duration * Fs / 2);
W_search_half_samples = round(W_search_duration * Fs / 2);
R_coarse_threshold = 0.3; % 成功粗略对齐的最小相关性


% --- 精细匹配与定位参数 ---
fine_match_yld_ch_idx = 1; % YLD 用于精细匹配的天线
fine_match_chj_ch_idx = 1; % CHJ 用于精细匹配的天线
findpeaks_noise_duration = 2e8; % 噪声的信号段持续时间
findpeaks_k = 5.0; % findpeaks 的阈值乘数（k * 噪声标准差）
findpeaks_min_dist_time = 1e-6; % 最小峰值距离
W_pulse_duration = 5e-6; % 精细相关性的脉冲窗口持续时间
W_fine_search_duration = 30e-6; % 精细搜索窗口在预测时间周围的持续时间
R_fine_threshold = 0.5; % 成功精细匹配的最小相关性
time_tolerance_s = 10e-6; % 允许的测量值与计算出的 TDOA 之间的最大允许差异

W_pulse_half_samples = round(W_pulse_duration * Fs / 2);
W_fine_search_half_samples = round(W_fine_search_duration * Fs / 2);
findpeaks_min_dist_samples = round(findpeaks_min_dist_time * Fs);
time_tolerance_samples = round(time_tolerance_s * Fs);

% --- 数据加载参数 ---
load_r_location_yld = 4.0e8; % YLD 数据读取起始位置
load_r_location_chj = 4.0e8; % CHJ 数据读取起始位置
load_signal_length = 1.5e8; % 要加载的信号长度

% --- 结果存储 ---
S_results = [];           % 用于存储 3D 位置
match_details = [];     % 用于存储详细的匹配信息（YLD 位置、CHJ 位置、相似度等）
%% Phase 1: 数据加载与预处理 =================================
% --- 加载数据 ---
fprintf('  正在加载 YLD 数据...\n');
yld_signals = cell(1, 3);
for ch = 1:3
    filename = fullfile(yld_data_path, sprintf('%sCH%d.dat', yld_base_filename, ch));
    fprintf('    正在读取 %s\n', filename);
    try
        yld_signals{ch} = read_signal(filename, load_signal_length, load_r_location_yld);
    catch ME
        error('无法读取 YLD 数据文件：%s\n%s', filename, ME.message);
    end
end

fprintf('  正在加载 CHJ 数据...\n');
chj_signals = cell(1, 3);
for ch = 1:3
    % 仔细构建 CHJ 文件名（处理空格）
    chj_fname_part1 = chj_base_filename(1:4); % "2024"
    chj_fname_part2 = chj_base_filename(6:8); % "822"
    chj_fname_part3 = chj_base_filename(10:end);% "85933.651462"
    filename = fullfile(chj_data_path, sprintf('%s %s %sCH%d.dat', chj_fname_part1, chj_fname_part2, chj_fname_part3, ch));
    fprintf('    正在读取 %s\n', filename);
    try
        % 重要：CHJ 数据加载相对于其自身的开始时间/位置
        chj_signals{ch} = read_signal(filename, load_signal_length, load_r_location_chj);
    catch ME
        error('无法读取 CHJ 数据文件：%s\n%s', filename, ME.message);
    end
end
fprintf('  数据加载完成。\n');

% --- 预处理（滤波） ---
fprintf('  正在应用滤波器...\n');
filtered_yld_signals = cell(1, 3);
filtered_chj_signals = cCHJ_FILE_OFFSETell(1, 3);
for ch = 1:3
    filtered_yld_signals{ch} = filter_xb(yld_signals{ch});
    filtered_chj_signals{ch} = filter_xb(chj_signals{ch});
end
clear yld_signals chj_signals; % 释放内存
fprintf('  滤波完成。\n');

%% Phase 2: 粗略对齐 =============================================
fprintf('Phase 2: 粗略对齐\n');

% --- 选择用于对齐的数据 ---
align_yld_signal = filtered_yld_signals{coarse_align_yld_ch_idx};
align_chj_signal = filtered_chj_signals{coarse_align_chj_ch_idx};

% --- 调整脉冲位置相对于加载的数据 ---
% 提供的脉冲位置需要相对于加载的信号段的起始点
anchor_loc_yld_relative = anchor_loc_yld_approx - load_r_location_yld + 1;
anchor_loc_chj_relative_approx = anchor_loc_chj_approx - load_r_location_chj + 1;

if anchor_loc_yld_relative <= 0 || anchor_loc_yld_relative > load_signal_length || ...
   anchor_loc_chj_relative_approx <= 0 || anchor_loc_chj_relative_approx > load_signal_length
    error('近似脉冲位置超出加载的数据段。调整加载参数或脉冲位置。');
end

% --- 定义 YLD 脉冲窗口 ---
yld_start = max(1, anchor_loc_yld_relative - W_coarse_half_samples);
yld_end = min(load_signal_length, anchor_loc_yld_relative + W_coarse_half_samples);
yld_anchor_window = align_yld_signal(yld_start : yld_end);
fprintf('  YLD 脉冲窗口定义（样本 %d 到 %d，相对）。\n', yld_start, yld_end);

% --- 定义 CHJ 搜索段 ---
chj_search_start = max(1, anchor_loc_chj_relative_approx - W_search_half_samples);
chj_search_end = min(load_signal_length, anchor_loc_chj_relative_approx + W_search_half_samples);
chj_search_segment = align_chj_signal(chj_search_start : chj_search_end);
fprintf('  CHJ 搜索段定义（样本 %d 到 %d，相对）。\n', chj_search_start, chj_search_end);

% --- 执行互相关 ---
fprintf('  正在执行粗略互相关（这可能需要一些时间）...\n');
[r_coarse, lags_coarse] = xcorr(chj_search_segment, yld_anchor_window, 'normalized');
fprintf('  相关性计算完成。\n');

% --- 寻找峰值并计算偏移量 ---
[max_r_coarse, idx_max_lag] = max(r_coarse);
best_lag_coarse = lags_coarse(idx_max_lag);

fprintf('  最大粗略相关性：%.4f\n', max_r_coarse);
if max_r_coarse < R_coarse_threshold
    warning('粗略对齐失败：最大相关性 (%.4f) 低于阈值 (%.4f)。', max_r_coarse, R_coarse_threshold);
    fprintf('  可能的问题:\n');
    fprintf('    - 脉冲选择不正确。\n');
    fprintf('    - CHJ 搜索窗口太小。\n');
    fprintf('    - 信号过于相似或噪声太大。\n');
    fprintf('    - 用于对齐的通道选择错误。\n');
    initial_offset_samples = input('粗略对齐失败。请输入手动偏移量猜测或按 Ctrl+C 终止：');
    if isempty(initial_offset_samples)
       error('粗略对齐失败且未提供手动偏移量。');
    end
else
    initial_offset_samples = best_lag_coarse;
    fprintf('  粗略对齐成功。\n');
    fprintf('  计算出的初始偏移量：%d 个样本。\n', initial_offset_samples);
    fprintf('  （这意味着 t_chj ≈ t_yld + %d 个样本）\n', initial_offset_samples);
end

% --- 验证图 ---
fprintf('  正在绘制粗略对齐验证图...\n');
figure('Name', '粗略对齐验证');
plot(load_r_location_chj:load_r_location_chj+load_signal_length-1, align_chj_signal, 'b');
hold on;
% 使用计算出的偏移量调整 YLD 时间轴
yld_time_axis_shifted = (load_r_location_yld:load_r_location_yld+load_signal_length-1) + initial_offset_samples;
plot(yld_time_axis_shifted, align_yld_signal, 'r');
xlabel('样本索引（对齐到 CHJ 时间基准）');
ylabel('滤波后振幅');
legend(sprintf('CHJ CH%d', coarse_align_chj_ch_idx), sprintf('YLD CH%d（偏移 %d)', coarse_align_yld_ch_idx, initial_offset_samples));
title(sprintf('粗略对齐验证（最大 R = %.3f)', max_r_coarse));
grid on;
% 在 CHJ 脉冲周围放大以目视检查对齐情况
xlim_center = anchor_loc_chj_relative_approx + load_r_location_chj -1;
xlim_range = W_search_half_samples / 10; % 放大因子
xlim([xlim_center - xlim_range, xlim_center + xlim_range]);
drawnow;
fprintf('  检查该图。如果对齐看起来合理，请按 Enter 继续。\n');
pause; % 等待用户确认

%% Phase 3: 脉冲识别 ==========================================
fprintf('Phase 3: 脉冲识别\n');

% --- 估计噪声水平 ---
noise_est_start_idx = 1; 
noise_est_end_idx = round(findpeaks_noise_duration * Fs);
noise_est_end_idx = min(noise_est_end_idx, load_signal_length);

noise_segment_yld = filtered_yld_signals{fine_match_yld_ch_idx}(noise_est_start_idx:noise_est_end_idx);
noise_segment_chj = filtered_chj_signals{fine_match_chj_ch_idx}(noise_est_start_idx:noise_est_end_idx);

% 使用中位数绝对偏差 (MAD) 以提高对异常值的鲁棒性
noise_std_yld = mad(noise_segment_yld, 1) * 1.4826;
noise_std_chj = mad(noise_segment_chj, 1) * 1.4826;

min_peak_height_yld = findpeaks_k * noise_std_yld;
min_peak_height_chj = findpeaks_k * noise_std_chj;

fprintf('  估计的噪声水平 (YLD CH%d): %.4f (阈值: %.4f)\n', fine_match_yld_ch_idx, noise_std_yld, min_peak_height_yld);
fprintf('  估计的噪声水平 (CHJ CH%d): %.4f (阈值: %.4f)\n', fine_match_chj_ch_idx, noise_std_chj, min_peak_height_chj);
% --- 寻找峰值 ---
fprintf('  正在寻找 YLD CH%d 信号中的峰值...\n', fine_match_yld_ch_idx);
[yld_peak_amps, yld_peak_locs_relative] = findpeaks(filtered_yld_signals{fine_match_yld_ch_idx}, ...
    'MinPeakHeight', min_peak_height_yld, ...
    'MinPeakDistance', findpeaks_min_dist_samples);
yld_peak_locs_absolute = yld_peak_locs_relative + load_r_location_yld - 1;  % 转换为绝对样本索引
fprintf('    在 YLD CH%d 中找到 %d 个峰值。\n', length(yld_peak_locs_absolute), fine_match_yld_ch_idx);

fprintf('  正在寻找 CHJ CH%d 信号中的峰值...\n', fine_match_chj_ch_idx);
[chj_peak_amps, chj_peak_locs_relative] = findpeaks(filtered_chj_signals{fine_match_chj_ch_idx}, ...
    'MinPeakHeight', min_peak_height_chj, ...
    'MinPeakDistance', findpeaks_min_dist_samples);
chj_peak_locs_absolute = chj_peak_locs_relative + load_r_location_chj - 1; % 转换为绝对样本索引
fprintf('    在 CHJ CH%d 中找到 %d 个峰值。\n', length(chj_peak_locs_absolute), fine_match_chj_ch_idx);

%% Phase 4: 精细匹配与定位 ================================
fprintf('Phase 4: 精细匹配与定位\n');

h_waitbar = waitbar(0, '正在处理精细匹配...');
num_yld_peaks = length(yld_peak_locs_absolute);

for i = 1:num_yld_peaks
    waitbar(i/num_yld_peaks, h_waitbar, sprintf('正在处理 YLD 峰值 %d/%d', i, num_yld_peaks));

    current_yld_loc_abs = yld_peak_locs_absolute(i);
    current_yld_loc_rel = current_yld_loc_abs - load_r_location_yld + 1; %  YLD 数据加载起始点相对位置

    % 1. 预测 CHJ 到达位置（绝对索引）
    predicted_chj_loc_abs = current_yld_loc_abs + initial_offset_samples + delta_t_theo_samples;

    % 2. 定义精细搜索窗口（绝对索引）
    chj_fine_search_start_abs = predicted_chj_loc_abs - W_fine_search_half_samples;
    chj_fine_search_end_abs = predicted_chj_loc_abs + W_fine_search_half_samples;

    % 3. 在绝对窗口内寻找候选 CHJ 峰值
    candidate_indices = find(chj_peak_locs_absolute >= chj_fine_search_start_abs & chj_peak_locs_absolute <= chj_fine_search_end_abs);

    if isempty(candidate_indices)
        continue; % 在搜索窗口内未找到候选者
    end

    candidate_chj_locs_abs = chj_peak_locs_absolute(candidate_indices);

    % 4. 对候选者执行精细匹配
    best_match_chj_loc_abs = -1;
    max_similarity = -Inf;
    best_t_gcc_samples = NaN; % 存储精细相关性的时间滞后

    % 提取 YLD 脉冲窗口（用于访问加载数据的相对索引）
    yld_pulse_start_rel = max(1, current_yld_loc_rel - W_pulse_half_samples);
    yld_pulse_end_rel = min(load_signal_length, current_yld_loc_rel + W_pulse_half_samples);
    yld_pulse_window = filtered_yld_signals{fine_match_yld_ch_idx}(yld_pulse_start_rel:yld_pulse_end_rel);
    % 应用窗口处理
    processed_yld_pulse = real(windowsignal(detrend(yld_pulse_window)));


    for j = 1:length(candidate_chj_locs_abs)
        current_chj_candidate_loc_abs = candidate_chj_locs_abs(j);
        current_chj_candidate_loc_rel = current_chj_candidate_loc_abs - load_r_location_chj + 1; % 相对索引

        % 检查相对索引是否在加载的 CHJ 数据范围内
        if current_chj_candidate_loc_rel <= 0 || current_chj_candidate_loc_rel > load_signal_length
            warning('Candidate CHJ peak %d is outside loaded data range. Skipping.', current_chj_candidate_loc_abs);
            continue;
        end

        % 提取 CHJ 脉冲窗口（相对索引）
        chj_pulse_start_rel = max(1, current_chj_candidate_loc_rel - W_pulse_half_samples);
        chj_pulse_end_rel = min(load_signal_length, current_chj_candidate_loc_rel + W_pulse_half_samples);
        chj_pulse_window = filtered_chj_signals{fine_match_chj_ch_idx}(chj_pulse_start_rel:chj_pulse_end_rel);
        % 应用窗口处理
        processed_chj_pulse = real(windowsignal(detrend(chj_pulse_window)));


        % 计算相似度（归一化互相关最大值）
        try
            % 确保用于相关性的长度相同
            len_yld = length(processed_yld_pulse);
            len_chj = length(processed_chj_pulse);
            len = min(len_yld, len_chj);

            if len < 10 % 任意的最小长度
                 similarity = -Inf;
                 t_gcc = NaN;
            else
                % 使用完整的 xcorr 寻找最大值和滞后
                [r_fine, lags_fine] = xcorr(processed_chj_pulse(1:len), processed_yld_pulse(1:len), 'normalized');
                [similarity, idx_max_fine] = max(r_fine);
                % 使用 cal_tau 来获得亚样本时间差
           
                 t_gcc = cal_tau(r_fine, lags_fine'); 
            end

            if similarity > max_similarity
                max_similarity = similarity;
                best_match_chj_loc_abs = current_chj_candidate_loc_abs;
                best_t_gcc_samples = t_gcc; % 存储这个最佳匹配的滞后
            end
        catch ME
            fprintf('警告: YLD 峰值 %d, CHJ 候选 %d 的精细相关性计算出错: %s\n', current_yld_loc_abs, current_chj_candidate_loc_abs, ME.message);
            continue;
        end
    end 

    % 5. 验证最佳匹配并进行定位
    if best_match_chj_loc_abs ~= -1 && max_similarity >= R_fine_threshold

        % ---  开始使用匹配对进行定位 ---
        % 测量出的时间差（样本数）
        delta_t_measured_samples = best_match_chj_loc_abs - current_yld_loc_abs;

        % 获取 CHJ 的 2D 角度，使用匹配的脉冲窗口
        % 需要提取所有 3 个 CHJ 通道在 best_match_chj_loc_abs 的信号
        best_match_chj_loc_rel = best_match_chj_loc_abs - load_r_location_chj + 1;
        chj_win_start_rel = max(1, best_match_chj_loc_rel - W_pulse_half_samples); 
        chj_win_end_rel = min(load_signal_length, best_match_chj_loc_rel + W_pulse_half_samples);

        % 检查窗口索引是否有效，然后提取
        if chj_win_start_rel > chj_win_end_rel || chj_win_start_rel > load_signal_length || chj_win_end_rel < 1
             fprintf('警告: CHJ 2D 角度计算中的窗口索引无效。跳过匹配。\n');
             continue;
        end

        processed_chj_win_ch1 = real(windowsignal(detrend(filtered_chj_signals{1}(chj_win_start_rel:chj_win_end_rel))));
        processed_chj_win_ch2 = real(windowsignal(detrend(filtered_chj_signals{2}(chj_win_start_rel:chj_win_end_rel))));
        processed_chj_win_ch3 = real(windowsignal(detrend(filtered_chj_signals{3}(chj_win_start_rel:chj_win_end_rel))));

        try
            % 计算 2D 定位
            [chj_azimuth, chj_elevation, chj_Rcorr, chj_t123] = get_2d_result_single_window(best_match_chj_loc_abs, processed_chj_win_ch1, processed_chj_win_ch2, processed_chj_win_ch3,'chj');
            if isempty(chj_azimuth) || chj_Rcorr < 0.3
                 fprintf('跳过匹配，因为 CHJ 2D 结果无效/较差 (R=%.2f)。\n', chj_Rcorr);
                 continue;
            end
        catch ME_2D
            fprintf('警告: 获取 CHJ 2D 结果时出错: %s。跳过匹配。\n', ME_2D.message);
            continue;
        end

        % --- 计算 YLD 角度从其脉冲 ---
         yld_win_start_rel = max(1, current_yld_loc_rel - W_pulse_half_samples);
         yld_win_end_rel = min(load_signal_length, current_yld_loc_rel + W_pulse_half_samples);
         processed_yld_win_ch1 = real(windowsignal(detrend(filtered_yld_signals{1}(yld_win_start_rel:yld_win_end_rel))));
         processed_yld_win_ch2 = real(windowsignal(detrend(filtered_yld_signals{2}(yld_win_start_rel:yld_win_end_rel))));
         processed_yld_win_ch3 = real(windowsignal(detrend(filtered_yld_signals{3}(yld_win_start_rel:yld_win_end_rel))));
         try
            [yld_azimuth, yld_elevation, yld_Rcorr, yld_t123] = get_2d_result_single_window(current_yld_loc_abs, processed_yld_win_ch1, processed_yld_win_ch2, processed_yld_win_ch3,'yld');
             if isempty(yld_azimuth) || yld_Rcorr < 0.3 % 添加有效性检查
                 fprintf('跳过匹配，因为 YLD 2D 结果无效/较差 (R=%.2f)。\n', yld_Rcorr);
                 continue;
             end
         catch ME_2D_YLD
            fprintf('警告: 获取 YLD 2D 结果时出错: %s。跳过匹配。\n', ME_2D_YLD.message);
            continue;
         end
         % --- 结束重新计算 YLD 角度 ---

        % 执行 3D 定位
        [R1_x, R1_y, R1_z] = sph2cart(deg2rad(90-yld_azimuth), deg2rad(yld_elevation),1); 
        [R2_x, R2_y, R2_z] = sph2cart(deg2rad(90-chj_azimuth), deg2rad(chj_elevation),1); 
        A1 = [R1_x, R1_y, R1_z];
        A2 = [R2_x, R2_y, R2_z];
        p_vec = chj_sit - yld_sit; % 从 YLD 到 CHJ 的向量

        C = cross(A1, A2);
        normC = norm(C);
        if normC < eps
            fprintf('向量 A1 和 A2 对于 YLD 峰值 %d 是平行的。跳过。\n', current_yld_loc_abs);
            continue; % 避免除以零或接近零
        end
        c_unit = C / normC;

        M = [A1', -A2', c_unit']; % 使用列向量形成矩阵

        if abs(det(M)) < 1e-6 % 检查矩阵是否奇异或接近奇异
             fprintf('矩阵 M 对于 YLD 峰值 %d 是奇异的或接近奇异的。跳过。\n', current_yld_loc_abs);
             continue;
        end

        % 求解 R1_value, R2_value, R3_value（距离/沿向量投影）
        try
            solution = M \ p_vec'; 
            R1_value = solution(1);
            R2_value = solution(2); % 注意：这里的 R2_value 对应于沿 -A2 的距离
            R3_value = solution(3); % 投影到 C 单位向量上
        catch ME_Solve
             fprintf('求解线性系统 M*sol=p 对于 YLD 峰值 %d 出错: %s。跳过。\n', current_yld_loc_abs, ME_Solve.message);
             continue;
        end

        % %          计算相对于 YLD 站点的 3D 源位置 S
        % %          原始公式看起来很复杂，让我们使用直接的向量定义：
        % %          源 S = YLD_sit + R1_value * A1
        % %          源 S = CHJ_sit + R2_value * A2  （错误，这里 R2_value 是沿 -A2 从 YLD 角度）
        % %          源 S = CHJ_sit - R2_value_actual * A2（其中 R2_value_actual 是从 CHJ 开始的距离）
        % %          如果可能，让我们使用交点公式，或者取平均值：
        %          S_est1 = yld_sit + R1_value * A1;
        %          S_est2 = chj_sit - (-R2_value) * A2; % R2_value from solve was distance along -A2
        % 
        %          % Check consistency or average? Original code had complex formulas based on R1/R2 values.
        %          % Let's use the first estimate for now, assuming A1, A2, p form a reasonable geometry
        %          sub_S = S_est1;
        %          % Or average if they are close:
        %          % if norm(S_est1 - S_est2) < threshold_dist
        %          %    sub_S = (S_est1 + S_est2) / 2;
        %          % else
        %          %    fprintf('Warning: Discrepancy between S_est1 and S_est2 for YLD peak %d. Using S_est1.\n', current_yld_loc_abs);
        %          %    sub_S = S_est1; % Or skip?
        %          % end
        R1_vec = R1_value * A1;
        R2_vec = (-R2_value) * A2;
        R3_vec = R3_value * c_unit;

        % --- !!! 重要警告 !!! ---
        % 对公式推导、变量定义（特别是R2_value在比较中的含义以及R2向量的来源）以及'+ p_vec'项的验证至关重要。

        if R1_value <= R2_value % - 检查R2_value的符号在这里是否重要！
            % 原始公式: sub_S = R1 + (R1_value / R2_value)*(R1_value / (R1_value + R2_value)) * R3;
            if R2_value == 0 || (R1_value + R2_value) == 0
                fprintf('警告: 原始公式1中的除以零已避免。跳过YLD峰值 %d。\n', current_yld_loc_abs);
                continue;
            end
            try
                weight_factor1 = (R1_value / R2_value) * (R1_value / (R1_value + R2_value));
                sub_S = yld_sit + R1_vec + weight_factor1 * R3_vec;
            catch ME_Formula1
                 fprintf('评估原始公式1时出错，跳过YLD峰值 %d: %s。\n', current_yld_loc_abs, ME_Formula1.message);
                 continue;
            end

        else % R1_value > R2_value
            % Use the second formula
            % 原始公式: sub_S = R2 - (R2_value / R1_value)* (R2_value / (R1_value + R2_value)) * R3 + p;
            if R1_value == 0 || (R1_value + R2_value) == 0
                 fprintf('警告: 原始公式2中的除以零已避免。跳过YLD峰值 %d。\n', current_yld_loc_abs);
                 continue;
            end
            try
                weight_factor2 = (R2_value / R1_value) * (R2_value / (R1_value + R2_value));
                % 严格按照公式应用，将R2解释为R2_vec（来自CHJ）和p解释为p_vec。
                % 如讨论所述，'R2 + p'可能旨在使结果相对于YLD原点。
                % 假设公式以某种方式直接产生绝对位置，或者
                % 假设公式产生相对于YLD原点的位置。添加yld_sit使其成为绝对位置。
                sub_S_relative_to_YLD_maybe = R2_vec - weight_factor2 * R3_vec + p_vec;
                sub_S = yld_sit + sub_S_relative_to_YLD_maybe;
            catch ME_Formula2
                 fprintf('评估原始公式2时出错，跳过YLD峰值 %d: %s。\n', current_yld_loc_abs, ME_Formula2.message);
                 continue;
            end
        end
        if ~isempty(sub_S) && all(isfinite(sub_S))
            % % 最终时间验证检查
            t_chj_prop = norm(sub_S - chj_sit) / c_mps;
            t_yld_prop = norm(sub_S - yld_sit) / c_mps;
            delta_t_prop_s = abs(t_chj_prop - t_yld_prop);

            % 使用基于峰值位置测量的时间差
            delta_t_measured_s = abs(delta_t_measured_samples / Fs);
            % 或者：如果认为通过精细相关延迟得到的时间差更准确，可以使用以下代码
            % delta_t_measured_s_corr = abs(best_t_gcc_samples / Fs); % Assuming t_gcc is in samples

            time_diff_error_s = abs(delta_t_prop_s - delta_t_measured_s);

            if time_diff_error_s <= time_tolerance_s
                fprintf('  匹配成功: YLD %d -> CHJ %d (R=%.3f). 位置: [%.1f, %.1f, %.1f]. 时间误差=%.2f微秒\n', ...
                        current_yld_loc_abs, best_match_chj_loc_abs, max_similarity, sub_S(1), sub_S(2), sub_S(3), time_diff_error_s*1e6);
                S_results = [S_results; sub_S];
                match_details = [match_details; struct(...
                    'yld_loc', current_yld_loc_abs, ...
                    'chj_loc', best_match_chj_loc_abs, ...
                    'similarity', max_similarity, ...
                    'delta_t_measured_s', delta_t_measured_s, ...
                    'delta_t_prop_s', delta_t_prop_s, ...
                    'loc_3d', sub_S)];
            else
                 fprintf('  匹配被拒绝（时间不一致）: YLD %d -> CHJ %d (R=%.3f). 时间误差=%.2f微秒 > %.2f微秒\n', ...
                         current_yld_loc_abs, best_match_chj_loc_abs, max_similarity, time_diff_error_s*1e6, time_tolerance_s*1e6);
            end
        end 
    end 
end 

close(h_waitbar);
fprintf('第4阶段: 精细匹配与定位完成。找到 %d 个有效的3D位置。\n', size(S_results, 1));

%% Phase 5: 结果存储与输出 =====================================
fprintf('Phase 5: 结果存储与输出\n');

% --- 保存结果 ---
results_filename = sprintf('localization_results_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
save(results_filename, 'S_results', 'match_details', 'yld_sit', 'chj_sit', 'initial_offset_samples');
fprintf('  结果已保存到 %s\n', results_filename);

% --- 绘制3D结果 ---
if ~isempty(S_results)
    figure('Name', '3D Localization Results');
    scatter3(S_results(:,1), S_results(:,2), S_results(:,3), 10, 'filled');
    hold on;
    scatter3(yld_sit(1), yld_sit(2), yld_sit(3), 100, 'r^', 'filled', 'DisplayName', 'YLD Station');
    scatter3(chj_sit(1), chj_sit(2), chj_sit(3), 100, 'bs', 'filled', 'DisplayName', 'CHJ Station');
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('3D Lightning Source Locations');
    legend show;
    axis equal; grid on;
    view(3); 
else
    fprintf('  No valid 3D locations found to plot.\n');
end

fprintf('Processing Finished.\n');