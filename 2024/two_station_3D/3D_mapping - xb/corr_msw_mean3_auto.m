%% 1. 主配置区域
clear;
clc;
close all;
c = 0.299792458;
fs = 200e6;
upsampling_factor = 50;
N = 3;
% --- 滤波器与阈值参数 ---
filter_band = [30e6, 80e6];
filter_order = 5;
noise_analysis_length = 1e8;
threshold_std_multiplier = 7;
% --- 分析参数 ---
window_lengths_to_run = [4096, 2048, 1024, 512, 256, 128, 64];
rcorr_thresholds_to_plot = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8];
base_output_dir = '2d_results—std—7——step_0.25';
% --- 信号段定义 ---
% 将每个数据集定义为 'segments' 结构体数组的一个元素。
% 脚本将自动遍历并处理这里定义的每一个数据段。

% 数据段 1: 引雷点
segments(1).name = 'yld';
segments(1).data_path = '..\\'; % .dat 数据文件的相对路径
segments(1).base_filename = '20240822165932.6610';
% --- 在下方数组中定义所有需要分析的起始位置 ---
segments(1).r_loctions = 3.65e8;
segments(1).signal_length = 0.5e7;
segments(1).d12 = 24.9586;
segments(1).d13 = 34.9335;
segments(1).d23 = 24.9675;
segments(1).angle12 = -110.8477;
segments(1).angle13 = -65.2405;
segments(1).angle23 = -19.6541;
segments(1).objective_type = 'yld'; % 'objective' 优化函数所需的类型标志

% % 数据段 2: 从化局
% segments(2).name = 'chj';
% segments(2).data_path = '..\\'; % .dat 数据文件的相对路径
% segments(2).base_filename = '2024 822 85933.651462';
% segments(2).r_loction = 3.95e8;
% segments(2).signal_length = 4e7;
% segments(2).d12 = 41.6496;
% segments(2).d13 = 36.9015;
% segments(2).d23 = 35.4481;
% segments(2).angle12 = -2.8381;
% segments(2).angle13 = 50.3964;
% segments(2).angle23 = 120.6568;
% segments(2).objective_type = 'chj'; % 'objective' 优化函数所需的类型标志

%% 2. 自动化处理
fprintf('===== 开始全自动雷电数据分析 =====\n');

% --- 最外层循环: 遍历每个定义的站点 ---
for s_idx = 1:length(segments)
    segment_template = segments(s_idx);

    % --- 遍历当前站点的每一个起始位置 (r_loction) ---
    for r_idx = 1:length(segment_template.r_loctions)
        current_segment = segment_template;
        current_segment.r_loction = segment_template.r_loctions(r_idx); % 设置当前起始位置

        % 生成能区分 r_loction 的文件夹名称
        segment_run_name = sprintf('%s_start_%.2fe8', current_segment.name, current_segment.r_loction / 1e8);

        fprintf('\n--- 正在处理: %s ---\n', segment_run_name);

        % --- 数据加载与预处理 (每个 r_loction 执行一次) ---
        fprintf('正在加载并滤波数据...\n');
        try
            ch1_file = fullfile(current_segment.data_path, [current_segment.base_filename, 'CH1.dat']);
            ch2_file = fullfile(current_segment.data_path, [current_segment.base_filename, 'CH2.dat']);
            ch3_file = fullfile(current_segment.data_path, [current_segment.base_filename, 'CH3.dat']);

            % 加载主信号
            ch1 = read_signal(ch1_file, current_segment.signal_length, current_segment.r_loction);
            ch2 = read_signal(ch2_file, current_segment.signal_length, current_segment.r_loction);
            ch3 = read_signal(ch3_file, current_segment.signal_length, current_segment.r_loction);

            % 滤波信号
            filtered_signal1 = filter_bp(ch1, filter_band(1), filter_band(2), filter_order);
            filtered_signal2 = filter_bp(ch2, filter_band(1), filter_band(2), filter_order);
            filtered_signal3 = filter_bp(ch3, filter_band(1), filter_band(2), filter_order);

            % 计算阈值
            noise = read_signal(ch1_file, noise_analysis_length, noise_analysis_length);
            filtered_noise = filter_bp(noise, filter_band(1), filter_band(2), filter_order);
            threshold = mean(filtered_noise) + threshold_std_multiplier * std(filtered_noise);
            fprintf('计算得到的噪声阈值: %f\n', threshold);

        catch ME
            fprintf(2, '为 %s 加载或滤波数据时出错。正在跳过...\n', segment_run_name);
            fprintf(2, '错误信息: %s\n', ME.message);
            continue;
        end

        % --- 内层循环: 遍历每个指定的窗口长度 ---
        for w_idx = 1:length(window_lengths_to_run)
            window_length = window_lengths_to_run(w_idx);

            % --- 寻找波峰  ---
            min_peak_dist = window_length/4;
            [all_peaks, all_locs] = findpeaks(filtered_signal1, 'MinPeakHeight', threshold, 'MinPeakDistance', min_peak_dist);
            num_peaks_total = numel(all_peaks);

            if isempty(all_peaks)
                fprintf('在 %s 未找到有效波峰。正在跳过...\n', segment_run_name);
                continue;
            end
            fprintf('找到 %d 个潜在波峰进行分析。\n', num_peaks_total);
            window = window_length * upsampling_factor;

            fprintf('\n> 正在使用窗口长度进行分析: %d\n', window_length);

            % --- 动态设置文件夹和文件名 ---
            output_dir = fullfile(base_output_dir, segment_run_name, ['window_', num2str(window_length)],['min_peak_dist_', num2str(min_peak_dist)]);
            if ~exist(output_dir, 'dir')
                mkdir(output_dir);
            end

            results_filename = sprintf('results_%s_win%d_threshold_%f.txt', ...
                current_segment.objective_type, window_length, threshold);
            results_filepath = fullfile(output_dir, results_filename);

            % --- 波峰处理 ---
            fileID = fopen(results_filepath, 'w');
            header = {'Start_loc', 'peak', 't12', 't13', 't23', 'cos_alpha_opt', ...
                'cos_beta_opt', 'Azimuth', 'Elevation', 'Rcorr', 't123'};
            fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', header{:});

            resultsTable = table('Size', [num_peaks_total, numel(header)], 'VariableTypes', ...
                repmat({'double'}, 1, numel(header)), 'VariableNames', header);
            results_count = 0;

            h_waitbar = waitbar(0, sprintf('数据: %s, 窗口: %d, 处理中...', ...
                replace(segment_run_name, '_', ' '), window_length));

            for pi = 1:num_peaks_total
                if mod(pi, 10) == 0
                    waitbar(pi / num_peaks_total, h_waitbar, sprintf('数据: %s, 窗口: %d, 正在处理波峰 %d/%d', ...
                        replace(segment_run_name, '_', ' '), window_length, pi, num_peaks_total));
                end

                idx = all_locs(pi);

                start_idx = idx - (window_length / 2 - 1);
                end_idx = idx + (window_length / 2);
                if start_idx <= 0 || end_idx > length(filtered_signal1)
                    continue;
                end

                signal1 = filtered_signal1(start_idx:end_idx);
                signal2 = filtered_signal2(start_idx:end_idx);
                signal3 = filtered_signal3(start_idx:end_idx);

                ch1_new = real(windowsignal(detrend(signal1)));
                ch2_new = real(windowsignal(detrend(signal2)));
                ch3_new = real(windowsignal(detrend(signal3)));

                [ch1_up, ch2_up, ch3_up] = deal(upsampling(ch1_new, upsampling_factor)', upsampling(ch2_new, upsampling_factor)', upsampling(ch3_new, upsampling_factor)');

                ch1_upsp = ch1_up(:, 2);
                ch2_upsp = ch2_up(:, 2);
                ch3_upsp = ch3_up(:, 2);

                [r12_gcc, lags12_gcc] = xcorr(ch1_upsp, ch2_upsp, 'normalized');
                [r13_gcc, lags13_gcc] = xcorr(ch1_upsp, ch3_upsp, 'normalized');
                [r23_gcc, lags23_gcc] = xcorr(ch2_upsp, ch3_upsp, 'normalized');

                R12_gcc = max(r12_gcc);
                R13_gcc = max(r13_gcc);
                R23_gcc = max(r23_gcc);

                t12_gcc = cal_tau(r12_gcc, lags12_gcc');
                t13_gcc = cal_tau(r13_gcc, lags13_gcc');
                t23_gcc = cal_tau(r23_gcc, lags23_gcc');

                t12 = t12_gcc * 0.1;
                t13 = t13_gcc * 0.1;
                t23 = t23_gcc * 0.1;

                cos_beta_0 = ((c*t13*current_segment.d12*sind(current_segment.angle12)) - (c*t12*sind(current_segment.angle13)*current_segment.d13)) / (current_segment.d13*current_segment.d12*sind(current_segment.angle12-current_segment.angle13));
                cos_alpha_0 = ((c*t12)/current_segment.d12 - cos_beta_0*cosd(current_segment.angle12))/sind(current_segment.angle12);

                if abs(cos_beta_0) > 1 || abs(cos_alpha_0) > 1, continue; end

                x0 = [cos_alpha_0, cos_beta_0];
                options = optimoptions('lsqnonlin', 'MaxIter', 1000, 'TolFun', 1e-6, 'Display', 'none');
                x = lsqnonlin(@(x) objective(x, t12, t13, t23, current_segment.objective_type), x0, [-1 -1], [1 1], options);

                cos_alpha_opt = x(1);
                cos_beta_opt = x(2);

                if abs(cos_alpha_opt) > 1 || abs(cos_beta_opt) > 1 || abs(cos_beta_opt/cos(atan2(cos_alpha_opt,cos_beta_opt))) > 1, continue; end

                Az = atan2(cos_alpha_opt, cos_beta_opt);
                El = acos(cos_beta_opt / cos(Az));
                Az_deg = rad2deg(Az);
                El_deg = rad2deg(El);
                if Az_deg < 0, Az_deg = Az_deg + 360; end

                t123 = t12 + t23 - t13;
                Rcorr = (R12_gcc + R13_gcc + R23_gcc) / 3;

                results_count = results_count + 1;
                current_result = {current_segment.r_loction + idx - window / 100, window / 100, t12, t13, t23, ...
                    cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr, t123};

                resultsTable(results_count, :) = current_result;
                fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', current_result{:});
            end

            close(h_waitbar);
            fclose(fileID);
            fprintf('结果已保存至 %s\n', results_filepath);

            resultsTable = resultsTable(1:results_count, :);

            % --- 绘图并保存图像 (现在会为每个Rcorr阈值循环)---
            if ~isempty(resultsTable)
                % 新增循环：遍历所有指定的Rcorr阈值
                for rc_idx = 1:length(rcorr_thresholds_to_plot)
                    current_rcorr_thresh = rcorr_thresholds_to_plot(rc_idx);
            
                    fprintf('>> 正在为 Rcorr > %.1f 生成图像...\n', current_rcorr_thresh);
            
                    % 使用当前的Rcorr阈值进行筛选
                    logicalIndex = abs(resultsTable.t123) < 1 & abs(resultsTable.Rcorr) > current_rcorr_thresh;
                    filteredTable = resultsTable(logicalIndex, :);
            
                    if ~isempty(filteredTable)
                        fig = figure('Visible', 'off'); % 在后台创建图像
                        Start_loc = filteredTable.Start_loc;
                        colorValues = (Start_loc - min(Start_loc)) / (max(Start_loc) - min(Start_loc));
                        if isnan(colorValues), colorValues = 0; end
            
                        scatter(filteredTable.Azimuth, filteredTable.Elevation, 1, colorValues, 'filled');
                        
                        % 更新标题以包含Rcorr信息
                        titleStr = sprintf('方位角 vs 俯仰角 (%s, 窗口: %d, Rcorr > %.1f)', ...
                            replace(segment_run_name, '_', ' '), window_length, current_rcorr_thresh);
                        title(titleStr);
                        
                        xlabel('方位角 (度)'); xlim([0, 360]); xticks(0:40:360);
                        ylabel('俯仰角 (度)'); ylim([0, 90]); yticks(0:10:90);
                        grid on;
                        cb = colorbar; ylabel(cb, '归一化起始位置');
                        colormap('jet'); caxis([0, 1]);
            
                        % 更新图像文件名以包含Rcorr信息
                        figure_filename = sprintf('Plot_Azimuth_Elevation_win%d_Rcorr_gt_%.1f.png', window_length, current_rcorr_thresh);
                        figure_filepath = fullfile(output_dir, figure_filename);
                        saveas(fig, figure_filepath);
                        close(fig); % 关闭图像以释放内存
                        fprintf('   图像已保存至: %s\n', figure_filepath);
                    else
                        fprintf('>> Rcorr > %.1f 条件下，经过筛选后没有可供绘制的数据点。\n', current_rcorr_thresh);
                    end
                end
            else
                fprintf('此窗口没有处理任何有效的波峰。\n');
            end
        end % --- 窗口长度循环结束 ---
    end % --- r_loction 循环结束 ---
end % --- 站点循环结束 ---

fprintf('\n===== 所有处理已完成。 =====\n');