%% Step1 读取引雷点的二维定位结果（需要条件筛选出合格的）
% 引入变量：位置，方位角，仰角
chj_signal_length = 5120;
match_signal_length = 6000;
yld_result_path = 'result_yld_window5120_3e8.txt';
start_signal_loc = 3.6e8;
mapping_start_signal_loc = 3.7e8;
end_signal_loc = 3.8e8;
step = 127200;
% 引入两个站的位置关系
yld_sit = [0, 0, 0];
chj_sit = [2003.7972, -7844.7836, -27];
% offsets = [-85000, -79000, -75000, -65000, -57000, -48000, -36000, -30000, -25000, -14000, -5000, 5000, 11000, 16000, 26000, 34000, 38000];
% yld相对于chj的位置
p = chj_sit-yld_sit;
dist = 8.0967e3; %单位：米
c = 0.299792458;
W = 30000; % 时间误差
offsets_init = -85438;
signal_length=step;
% 所有信号的开始位置
all_start_signal_loc = start_signal_loc:step:end_signal_loc;
all_S_results = [];   % 存储DTOA优化后的结果
all_match_results = []; % 存储DTOA优化后的匹配信息
% YLD站内天线局部坐标
yld_ant_local_coords = [
    0,   0,   0;  % 天线1
    -23.3246,  -8.8824,   0;  % 天线2
    -31.7222,  14.6305,   0   % 天线3
    ];
% CHJ站内天线局部坐标
chj_ant_local_coords = [
    0,   0,   0;  % 天线1
    -2.0622, 41.5985,   0;  % 天线2
    28.4316, 23.5237,   0   % 天线3
    ];

% 将局部天线坐标转换为全局坐标
yld_ant_global_coords = yld_ant_local_coords + yld_sit;
chj_ant_global_coords = chj_ant_local_coords + chj_sit;

baseline_ant_indices = [
    1, 2; % 基线12
    1, 3; % 基线13
    2, 3  % 基线23
    ];

% --- 所有基线的天线全局坐标 (6条站内 + 1条站间) ---
all_baseline_definitions = zeros(7, 6);
all_baseline_definitions(1, :) = [yld_ant_global_coords(1,:), yld_ant_global_coords(2,:)];
all_baseline_definitions(2, :) = [yld_ant_global_coords(1,:), yld_ant_global_coords(3,:)];
all_baseline_definitions(3, :) = [yld_ant_global_coords(2,:), yld_ant_global_coords(3,:)];
all_baseline_definitions(4, :) = [chj_ant_global_coords(1,:), chj_ant_global_coords(2,:)];
all_baseline_definitions(5, :) = [chj_ant_global_coords(1,:), chj_ant_global_coords(3,:)];
all_baseline_definitions(6, :) = [chj_ant_global_coords(2,:), chj_ant_global_coords(3,:)];
% 站间基线 (YLD参考天线 - CHJ参考天线)
all_baseline_definitions(7, :) = [yld_sit, chj_sit];
all_residuals_history = []; 
% --- DTOA测量不确定度  ---
sigmas_ns = [
    4.22; 5.54; 7.95; % YLD 站内基线 DTOA 不确定度 (ns)
    3.61; 13.89; 39.28; % CHJ 站内基线 DTOA 不确定度 (ns)
    3000         % 站间基线 DTOA 不确定度 (ns)
    ];



for j = 1:numel(all_start_signal_loc)-1
    start_read_loc_yld = all_start_signal_loc(j);
    end_read_loc_yld = all_start_signal_loc(j+1);
    if  start_read_loc_yld < mapping_start_signal_loc || end_read_loc_yld > end_signal_loc
        continue
    end
    % 记录处理的位置
    fprintf('正在处理的信号位置：%d -- %d \n', start_read_loc_yld, end_read_loc_yld);
    [yld_start_loc, yld_t12, yld_t13,yld_t23,yld_azimuth, yld_elevation, yld_Rcorr, yld_t123] = read_result(yld_result_path,start_read_loc_yld, end_read_loc_yld);
    yld_ch1 =read_signal('..\\20240822165932.6610CH1.dat',signal_length,start_read_loc_yld);
    chj_ch1 =read_signal('..\\2024 822 85933.651462CH1.dat',signal_length,start_read_loc_yld+ 34151156 - offsets_init-(j-1)*100);
    chj_ch2 =read_signal('..\\2024 822 85933.651462CH2.dat',signal_length,start_read_loc_yld+ 34151156 - offsets_init-(j-1)*100);
    chj_ch3 =read_signal('..\\2024 822 85933.651462CH3.dat',signal_length,start_read_loc_yld+ 34151156 - offsets_init-(j-1)*100 +215/5);
    filtered_yld_signal1 = filter_bp(yld_ch1,30e6,80e6,5);
    filtered_chj_signal1 = filter_bp(chj_ch1,30e6,80e6,5);
    filtered_chj_signal2 = filter_bp(chj_ch2,30e6,80e6,5);
    filtered_chj_signal3 = filter_bp(chj_ch3,30e6,80e6,5);
    S_results = [];
    match_results = struct('yld_start_loc', {}, 'chj_loc', {}, 'r_gccs', {}, 'chj_azimuth',{},'chj_elevation',{},'dlta',{});
    h = waitbar(0, 'Processing...');
    %% Step2 根据引雷点的信号窗口得到匹配到的从化局的信号
    for i =1 :numel(yld_start_loc)
        sub_S_results = [];
        sub_R_gccs = [];
        dltas = [];
        sub_chj_locs = [];
        sub_chj_azimuth = [];
        sub_chj_elevation = [];
        sub_chj_dtoa_t = [];
        yld_chj_dlta_ts = [];
        yld_chj_dlta_Ts = [];
        waitbar(i/numel(yld_start_loc), h, sprintf('位置：%d -- %d ；进度： %.2f%%', start_read_loc_yld, end_read_loc_yld, i/numel(yld_start_loc)*100));
        if yld_Rcorr(i) < 0.3 && yld_t123(i) > 1
            continue
        end
        % 转换绝对位置到相对位置
        yld_signal_start_loc = yld_start_loc(i) - start_read_loc_yld;
        start_read_loc_chj = yld_signal_start_loc;
        if yld_signal_start_loc+chj_signal_length > signal_length || yld_signal_start_loc + 1 < 0
            continue
        end
        processed_yld_signal = filtered_yld_signal1(yld_signal_start_loc+1 : yld_signal_start_loc+chj_signal_length);
        processed_yld_signal = real(windowsignal(detrend(processed_yld_signal)));

        % 读取从化局与引雷点同位置前后加6000点的信号
        if start_read_loc_chj-match_signal_length+1 <= 0 || start_read_loc_chj+match_signal_length+chj_signal_length > signal_length
            continue;
        end
        chj_match_signal1 = filtered_chj_signal1(start_read_loc_chj-match_signal_length+1:start_read_loc_chj+match_signal_length+chj_signal_length);
        chj_match_signal2 = filtered_chj_signal2(start_read_loc_chj-match_signal_length+1:start_read_loc_chj+match_signal_length+chj_signal_length);
        chj_match_signal3 = filtered_chj_signal3(start_read_loc_chj-match_signal_length+1:start_read_loc_chj+match_signal_length+chj_signal_length);

        % 寻找峰值
        [peaks, locs] = findpeaks(chj_match_signal1, 'MinPeakHeight', 5, 'MinPeakDistance', 1024);
        all_locs = locs;
        % 遍历所有峰值
        num_peaks = numel(all_locs);
        if num_peaks == 0
            continue;
        end

        for pi = 1:num_peaks
            idx = all_locs(pi);
            if idx - (chj_signal_length / 2 - 1) <= 0 || idx + (chj_signal_length / 2) > match_signal_length*2+chj_signal_length
                continue;
            end
            processed_chj_signal1 = chj_match_signal1(idx - (chj_signal_length / 2)+ 1:idx + (chj_signal_length / 2));
            processed_chj_signal2 = chj_match_signal2(idx - (chj_signal_length / 2)+ 1:idx + (chj_signal_length / 2));
            processed_chj_signal3 = chj_match_signal3(idx - (chj_signal_length / 2)+ 1:idx + (chj_signal_length / 2));

            processed_chj_signal1 = real(windowsignal(detrend(processed_chj_signal1)));
            processed_chj_signal2 = real(windowsignal(detrend(processed_chj_signal2)));
            processed_chj_signal3 = real(windowsignal(detrend(processed_chj_signal3)));

            [r_gcc, lags_gcc] = xcorr(processed_yld_signal, processed_chj_signal1, 'normalized');
            R_gcc = max(r_gcc);
            t_gcc = cal_tau(r_gcc, lags_gcc');
            chj_start_idx = idx - (chj_signal_length / 2)+ 1 - t_gcc;
            chj_end_idx  = idx + (chj_signal_length / 2) - t_gcc;
            if chj_end_idx > length(chj_match_signal1) || chj_end_idx < 1 || chj_start_idx < 1 || chj_start_idx > length(chj_match_signal1)
                continue
            end
            processed_chj_signal1 = chj_match_signal1(chj_start_idx:chj_end_idx);
            processed_chj_signal2 = chj_match_signal2(chj_start_idx:chj_end_idx);
            processed_chj_signal3 = chj_match_signal3(chj_start_idx:chj_end_idx);

            processed_chj_signal1 = real(windowsignal(detrend(processed_chj_signal1)));
            processed_chj_signal2 = real(windowsignal(detrend(processed_chj_signal2)));
            processed_chj_signal3 = real(windowsignal(detrend(processed_chj_signal3)));
            [chj_start_loc, chj_t12, chj_t13, chj_t23, chj_azimuth, chj_elevation, chj_Rcorr, chj_t123] = get_2d_result_single_window(start_read_loc_chj,processed_chj_signal1,processed_chj_signal2,processed_chj_signal3,'chj');
            if chj_start_loc == 0 && chj_Rcorr < 0.3 && chj_t123 > 1
                continue
            end
            [R1_x, R1_y, R1_z] = sph2cart(deg2rad(90-yld_azimuth(i)), deg2rad(yld_elevation(i)),1);
            [R2_x, R2_y, R2_z] = sph2cart(deg2rad(90-chj_azimuth), deg2rad(chj_elevation),1);
            A1 = [R1_x, R1_y, R1_z];
            A2 = [R2_x, R2_y, R2_z];
            C = cross(A1, A2);
            if norm(C) < eps
                continue;  % 避免除以零
            end
            c_unit = C  / norm(C);  % 单位向量
            M = [A1(1), -A2(1), c_unit(1);
                A1(2), -A2(2), c_unit(2);
                A1(3), -A2(3), c_unit(3)];
            % 使用克莱姆法则求R1,R2,R3的标量
            detM = det(M);
            detR1 = det([p', M(:,2), M(:,3)]);
            detR2 = det([M(:,1), p', M(:,3)]);
            detR3 = det([M(:,1), M(:,2), p']);
            R1_value = detR1 / detM;
            R2_value = detR2 / detM;
            R3_value = detR3 / detM;
            R1 = R1_value * A1;
            R2 = R2_value * A2;
            R3 = R3_value/norm(C)* C;
            if R1_value <= R2_value
                % 使用第一个公式
                sub_S = R1 + (R1_value / R2_value)*(R1_value / (R1_value + R2_value)) * R3;
            else
                % 使用第二个公式
                sub_S = R2 - (R2_value / R1_value)* (R2_value / (R1_value + R2_value))  * R3 + p;
            end
            if sub_S(3) < 0
                sub_S = -sub_S;
            end
            if ~isempty(sub_S)
                t_chj = sqrt(sum((sub_S - chj_sit).^2))/c;
                t_yld = sqrt(sum((sub_S - yld_sit).^2))/c;
                dlta_t = t_yld-t_chj;
                dlta_T = (match_signal_length - chj_start_idx)*5;
                dlta = abs(dlta_t-dlta_T);
                if dlta <= W
                    dltas = [dltas;dlta];
                    yld_chj_dlta_ts = [yld_chj_dlta_ts;dlta_t];
                    yld_chj_dlta_Ts = [yld_chj_dlta_Ts;dlta_T];
                    sub_chj_dtoa_t = [sub_chj_dtoa_t;[chj_t12,chj_t13,chj_t23]];
                    sub_S_results = [sub_S_results; sub_S];
                    sub_R_gccs = [sub_R_gccs;R_gcc];
                    sub_chj_locs = [sub_chj_locs;chj_start_idx];
                    sub_chj_azimuth = [sub_chj_azimuth;chj_azimuth];
                    sub_chj_elevation = [sub_chj_elevation;chj_elevation];
                end
            end
        end
        [max_R_gcc, max_R_gcc_index] = max(sub_R_gccs);
        if isempty(max_R_gcc)
            continue
        end

        %% Step 5: 差分到达时间 (DTOA) 技术
        % --- 收集所有测量的DTOA值 (单位: ns) ---
        all_measured_dtoas_ns = zeros(7,1);

        all_measured_dtoas_ns(1:3) = [yld_t12(i),yld_t13(i),yld_t23(i)];  % YLD站内 (ns)
        all_measured_dtoas_ns(4:6) = sub_chj_dtoa_t(max_R_gcc_index);% CHJ站内 (ns)
        % 站间DTOA (ns)
        all_measured_dtoas_ns(7) = yld_chj_dlta_Ts(max_R_gcc_index);

        % --- DTOA 优化 ---
        S_initial = sub_S_results(max_R_gcc_index,:); % 使用三角测量结果作为初值

        objective_fun = @(S_vec) dtoa_objective_function_ns(S_vec, ...
            all_baseline_definitions, ...
            all_measured_dtoas_ns, ...
            sigmas_ns, ...
            c);

        options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', ...
            'Display', 'off', 'TolFun', 1e-3, 'TolX', 1e-3, ...
            'MaxIterations', 50, 'StepTolerance', 1e-3);
        S_optimized = S_initial; % 默认值
        optimization_successful = false;
        try
            [S_optimized_temp, ~, ~, exitflag, ~] = lsqnonlin(objective_fun, S_initial, [], [], options);
            if exitflag > 0 % 检查优化是否成功收敛
                S_optimized = S_optimized_temp;
                optimization_successful = true;
            else
                fprintf('优化未收敛 (exitflag=%d) for YLD event %d (sample loc %d)\n', exitflag, i, current_yld_start_abs_samples);
            end
        catch ME_optim
            fprintf('DTOA优化错误 for YLD event %d (sample loc %d): %s\n', i, yld_start_loc(i), ME_optim.message);
        end

        if optimization_successful
            % 用最终的优化解 S_optimized 再次调用目标函数，获取最终的理论DTOA
            final_theoretical_dtoas_ns = dtoa_objective_function_ns(S_optimized, ...
                all_baseline_definitions, ...
                all_measured_dtoas_ns, ...
                sigmas_ns, ...
                c);
            % 残差 = 理论DTOA - 测量DTOA
            final_residuals_ns = final_theoretical_dtoas_ns;

            all_residuals_history = [all_residuals_history; final_residuals_ns(:)'];

        end

        % --- 使用优化后的S重新进行时间校验 ---
        % 理论传播时间 (ns)
        t_chj_optimized_ns = norm(S_optimized - chj_sit) / c;
        t_yld_optimized_ns = norm(S_optimized - yld_sit) / c;
        % 理论站间DTOA (ns),
        theoretical_inter_dtoa_optimized_ns = t_yld_optimized_ns - t_chj_optimized_ns;

        % 最终校验差 (ns)
        final_dlta_check_ns = abs(theoretical_inter_dtoa_optimized_ns - yld_chj_dlta_Ts(max_R_gcc_index));

        if final_dlta_check_ns <= W && optimization_successful
            if S_optimized(3) < 0
                S_optimized = -S_optimized;
            end

            all_S_results = [all_S_results; S_optimized];
            match_info_dtoa = struct(...
                'yld_start_loc', yld_start_loc(i), ...
                'chj_loc', sub_chj_locs(max_R_gcc_index)+ start_read_loc_yld + 34151156 - offsets_init-(j-1)*100 +start_read_loc_chj-match_signal_length+1, ...
                'chj_azimuth', sub_chj_azimuth(max_R_gcc_index), ...
                'chj_elevation', sub_chj_elevation(max_R_gcc_index), ...
                'r_gccs', max_R_gcc, ...
                'dlta', final_dlta_check_ns, ...
                'S_initial_triangulation', S_initial ...
                );
            all_match_results = [all_match_results; match_info_dtoa];
        else
            % 如果DTOA优化失败或校验不通过，保存三角测量结果
            if dltas(max_R_gcc_index) <= W
                if S_initial(3) < 0
                    S_initial = -S_initial;
                end
                all_S_results = [all_S_results; S_initial];
                match_info_dtoa = struct(...
                    'yld_start_loc', yld_start_loc(i), ...
                    'chj_loc', sub_chj_locs(max_R_gcc_index)+ start_read_loc_yld + 34151156 - offsets_init-(j-1)*100 +start_read_loc_chj-match_signal_length+1, ...
                    'chj_azimuth', sub_chj_azimuth(max_R_gcc_index), ...
                    'chj_elevation', sub_chj_elevation(max_R_gcc_index), ...
                    'r_gccs', max_R_gcc, ...
                    'dlta', dltas(max_R_gcc_index), ...
                    'S_initial_triangulation', S_initial ...
                    );
                all_match_results = [all_match_results; match_info_dtoa];
            end
        end
    end
    close(h);
end



if ~isempty(all_residuals_history)

    baseline_names = {
        'YLD 1-2'; 'YLD 1-3'; 'YLD 2-3'; ...
        'CHJ 1-2'; 'CHJ 1-3'; 'CHJ 2-3'; ...
        'Inter-Station YLD-CHJ'
    };
    
    num_baselines = size(all_residuals_history, 2);
    num_events = size(all_residuals_history, 1);
    
    fprintf('\n--- 逐条基线不确定度验证 ---\n');
    
    % 创建一个新的 figure 用于绘制所有基线的残差直方图
    figure;
    sgtitle(sprintf('所有基线的残差分布 (基于 %d 个事件)', num_events)); % 为整个图添加总标题
    
    for k_base = 1:num_baselines
        % 提取当前基线的所有残差
        current_residuals = all_residuals_history(:, k_base);
        
        % 计算残差的均值和标准差
        mean_residual = mean(current_residuals);
        std_dev_residual = std(current_residuals);
        
        % 获取预设的不确定度值
        assumed_sigma = sigmas_ns(k_base);
        
        % --- 在命令窗口打印数值结果 ---
        fprintf('\n--- 基线 #%d: %s ---\n', k_base, baseline_names{k_base});
        fprintf('预设的不确定度 (σ): %.2f ns\n', assumed_sigma);
        fprintf('从 %d 个事件的残差中计算出的统计结果:\n', num_events);
        fprintf('  - 残差均值: %.2f ns\n', mean_residual);
        fprintf('  - 残差标准差: %.2f ns\n', std_dev_residual);
        
        % --- 在子图中绘制直方图和正态分布拟合曲线 ---
        subplot(3, 3, k_base); % 创建一个 3x3 的子图布局，并激活第 k_base 个
        histogram(current_residuals, 50, 'Normalization', 'pdf', 'EdgeColor', 'none');
        hold on;
        
        % 拟合一个正态分布曲线进行对比
        x_range = linspace(min(current_residuals), max(current_residuals), 100);
        pdf_fit = normpdf(x_range, mean_residual, std_dev_residual);
        plot(x_range, pdf_fit, 'r-', 'LineWidth', 1.5);
        
        title(baseline_names{k_base});
        xlabel('残差 (ns)');
        ylabel('概率密度');
        grid on;
        legend('残差分布', '正态拟合');
        hold off;
    end
end
% plot_3d;

% --- DTOA 目标函数 ---
function errors = dtoa_objective_function_ns(S_guess_xyz, all_baseline_ant_coords, all_measured_dtoas_ns, all_sigmas_ns, c_mpns)
num_baselines = size(all_baseline_ant_coords, 1);
errors = zeros(num_baselines, 1);
S_guess_xyz = S_guess_xyz(:)'; % 确保是行向量 [1x3]

for k_base = 1:num_baselines
    ant1_pos = all_baseline_ant_coords(k_base, 1:3);
    ant2_pos = all_baseline_ant_coords(k_base, 4:6);

    dist_to_ant1 = norm(S_guess_xyz - ant1_pos);
    dist_to_ant2 = norm(S_guess_xyz - ant2_pos);

    % 理论 DTOA (ns)
    theoretical_dtoa_ns = (dist_to_ant1 - dist_to_ant2) / c_mpns;

    errors(k_base) = (theoretical_dtoa_ns - all_measured_dtoas_ns(k_base)) / all_sigmas_ns(k_base);
end
end
