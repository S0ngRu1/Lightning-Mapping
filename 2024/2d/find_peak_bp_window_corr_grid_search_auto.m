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




function tau = cal_tau(R, lag)
    % 从数据中找到y的最大值及其索引
    [~, max_index] = max(R);
    tau = lag(max_index,1);
end

% function tau = cal_tau(R, lag) % 新的、高精度的版本
%     [~, max_idx] = max(R);
% 
%     % 确保峰值不在数组的边缘，否则无法取到3个点
%     if max_idx == 1 || max_idx == length(R)
%         tau = lag(max_idx);
%         return;
%     end
% 
%     % 提取峰值点 (y2) 和它左右相邻的两个点 (y1, y3)
%     y1 = R(max_idx - 1);
%     y2 = R(max_idx);
%     y3 = R(max_idx + 1);
% 
%     % 抛物线顶点横坐标的偏移量公式： p = (y1 - y3) / (2 * (y1 - 2*y2 + y3))
%     % p 是相对于中心点 max_idx 的亚采样偏移量
%     % 注意：要处理分母为0或非常小的情况，避免计算错误
%     denominator = 2 * (y1 - 2*y2 + y3);
%     if abs(denominator) < 1e-9
%         p = 0; % 如果分母太小（例如，平顶），则不进行偏移
%     else
%         p = (y1 - y3) / denominator;
%     end
%     
%     % 计算最终的精确时延
%     % lag是等差数列，可以直接用 p 乘以步长
%     time_step = lag(2) - lag(1);
%     tau = lag(max_idx) + p * time_step;
% end

function fitted_peak_x = fitpeak(data,peak_index)
if peak_index+10 < 10240 && peak_index-10 > 0
    fit_range = (peak_index + (-10:10))';
elseif peak_index+6 < 10240 && peak_index-6 > 0
    fit_range = (peak_index + (-6:6))';
elseif peak_index+2 < 10240 && peak_index-2 > 0
    fit_range = (peak_index + (-2:2))';
else
    fitted_peak_x = peak_index;
    return;
end
fit_values = data(fit_range);
coefficients = polyfit(fit_range, fit_values, 2);
fit_indices_curve = linspace(min(fit_range), max(fit_range), 1000);
fit_values_curve = polyval(coefficients, fit_indices_curve);
% 绘制原始数据和拟合曲线
% figure;
% plot(1:length(data),data)
% plot(1:length(data),data, 'b', fit_indices_curve, fit_values_curve, 'r--');
% legend('原始数据', '拟合曲线');
% xlabel('y的索引');
% ylabel('y的值');
[~, max_index_fit] = max(fit_values_curve);
fitted_peak_x = fit_indices_curve(1,max_index_fit);
end




function delta_t = delta_t(tij,tij_obs)
    delta_t = tij - tij_obs;
end



% 定义计算τij的理想值τ_ij^obs的函数
function tau_ij_obs = calculate_tau_obs(cos_alpha, cos_beta, type)
    % 初始化输出变量
    tau_ij_obs = zeros(1, 3);

    % 根据 type 参数选择不同的参数集
    if strcmp(type, 'chj') % 从化局
        angle12 = -2.8381;
        angle13 = 28.2006;
        angle23 = 87.3358;
        d12 = 41.6496;
        d13 = 48.5209;
        d23 = 25.0182;
    elseif strcmp(type, 'yld') % 引雷场
        angle12 = -110.8477;
        angle13 = -65.2405;
        angle23 = -19.6541;
        d12 = 24.9586;
        d13 = 34.9335;
        d23 = 24.9675;
    else
        error('未知的类型：%s', type);
    end

    % 使用式(3)计算τij的理想值τ_ij^obs
    tau_ij_obs(1) = (cos_alpha * sind(angle12) + cos_beta * cosd(angle12)) * d12 / 0.299792458;
    tau_ij_obs(2) = (cos_alpha * sind(angle13) + cos_beta * cosd(angle13)) * d13 / 0.299792458;
    tau_ij_obs(3) = (cos_alpha * sind(angle23) + cos_beta * cosd(angle23)) * d23 / 0.299792458;
end



%% 设计巴特沃斯带通滤波器
function filtered_signal = filter_bp(signal,f1,f2,order)
    Fs = 200e6;
    fn = Fs/2;
    Wn = [f1 f2]/fn;
    [b,a] = butter(order,Wn); 
    filtered_signal = filtfilt(b,a,signal);

end



function filtered_signal = filter_fft(sig,f1,f2)
    y=fft(sig);%傅里叶变换得到一个复数
    fs = 200e6;
    n = length(y);
    %创建一个长度与输入信号 y 相同的零向量 yy；
    yy=zeros(1,length(y));
    % 使用 for 循环遍历信号 y 的每个采样点（m 表示当前的采样点索引，从0到 N-1）；
    for m=1:n-1
    %     判断当前采样点对应的频率是否在 8Hz 到 15Hz 范围内，如果在该范围内，则将对应的 yy 值置为0，表示该频率的信号被滤除；
        if m*(fs/n)<f1 || m*(fs/n)>f2 %将奈奎斯特之后的频率也滤除点掉
            yy(m+1)=0;
        else
    %         如果当前采样点对应的频率不在 8Hz 到 15Hz 范围内，则将 yy 的值保持为原始信号 y 的值。
            yy(m+1)=y(m+1);
        end
    end %将频率为8Hz-15Hz的信号的幅值置0
    filtered_signal=ifft(yy)';
    
end


% 设计巴特沃斯带通滤波器
function filtered_signal = filtersignal(signal,f1,f2,order,fs)
     % 滤波器通带下边界频率f1 滤波器通带上边界频率f2  滤波器阶数order
     % 创建滤波器对象
     filter = designfilt('bandpassiir', 'FilterOrder', order, 'HalfPowerFrequency1', f1, 'HalfPowerFrequency2', f2, 'SampleRate', fs);
     filtered_signal = filtfilt(filter,signal);
end





%函数：遍历整个信号，找到微尺度窗口下相关系数大于0.8的窗口
function correlated_windows = find_correlated_windows(signal1, signal2, signal3, window_size, threshold, Fs, N)
    % 窗口数量
    num_windows = N - window_size + 1;  
    % 存储相关系数大于阈值的窗口
    correlated_windows = [];
    for i = 1:num_windows
        % 提取当前窗口的数据
        window1 = signal1(i:i+window_size-1);
        window2 = signal2(i:i+window_size-1);
        window3 = signal3(i:i+window_size-1);
        [~,R12,~] = gccphat(window1,window2, Fs);
        [~,R13,~] = gccphat(window1,window3, Fs);
        [~,R23,~] = gccphat(window2,window3, Fs);
        max_R12 = maxvalue(R12);
        max_R13 = maxvalue(R13);
        max_R23 = maxvalue(R23);
        
        % 如果相关系数大于阈值，将窗口添加到结果列表中
        if max_R12 > threshold && max_R13 > threshold && max_R23 > threshold
            correlated_windows = [correlated_windows; i];
        end
    end
end

function filtered_signal = rfi_filter(ori_signal,sub_signal_length)
    filtered_signal = [];
    subsignal_starts = 1:sub_signal_length/2:length(ori_signal);
    for i = 1:length(subsignal_starts)
        if subsignal_starts(i) + sub_signal_length - 1 > length(ori_signal)
            continue
        end
        subsignal = ori_signal(subsignal_starts(i):subsignal_starts(i)+sub_signal_length-1);
        windowed_signal = window_plus(sub_signal_length,subsignal);
        sub_filtered_signal = datafilter(windowed_signal);
        filtered_signal = [filtered_signal; sub_filtered_signal(sub_signal_length*0.25+1:sub_signal_length*0.75)];
    end
end


%函数：寻找信号的最大峰值
function peaks = find_max_peaks(signal,threshold)
    % 找到信号中的峰值
    [pks,locs] = findpeaks(signal);
    % 根据阈值筛选峰值
    selectedLocs = locs(pks > threshold);
    % 获取过滤后的每个峰值的x值
    
    peaks =selectedLocs;
end


%函数：寻找信号的峰值
function peaks = find_peaks(signal,threshold)
    % 找到信号中的峰值
    [pks,locs] = findpeaks(signal);
    % 根据阈值筛选峰值
    peaks = locs(pks > threshold);
end


function matched_peaks_x = match_peaks(peaks1,peaks2,peaks3)
    matched_peaks_x = []; % 存储匹配峰值的x值矩阵
    for i = 1:numel(peaks1)
        curr_peak1 = peaks1(i);
        % 检查peaks2和peaks3中是否存在与peaks1对应的峰值且x值的差不大于4
        idx_peak2 = find(abs(peaks2 - curr_peak1) <= 10);  % 获取peaks2中匹配峰值的索引
        idx_peak3 = find(abs(peaks3 - curr_peak1) <= 10);  % 获取peaks3中匹配峰值的索引
        % 检查是否找到了匹配的峰值
        if ~isempty(idx_peak2) && ~isempty(idx_peak3)
            matched_peaks_x = [matched_peaks_x; [curr_peak1, peaks2(idx_peak2(1)), peaks3(idx_peak3(1))]];% 添加匹配峰值的x值矩阵
        end
    end
end


function max_index = maxindex(vector)
    % 提取实部部分
    
    max_value = max(vector);
    % 找到最大值对应的索引
    max_index = find(vector == max_value);
end


function mswed_signal = msw_signal(signal , peak_x ,length)
      % 找到峰值的 x 值在信号中的索引
    left_idx = max(peak_x - length+1, 1);  % 确定左边界的索引
    right_idx = min(peak_x + length, 10240);  % 确定右边界的索引
    mswed_signal = signal(left_idx:right_idx);  % 提取以中心 x 值为中心的左右40个采样点

end

% 定义目标函数
function F = objective(x)
    % 提取待优化的变量
    cos_alpha = x(1);
    cos_beta = x(2);

    % 计算τij的理想值τ_ij^obs
    tau_ij_obs = calculate_tau_obs(cos_alpha, cos_beta);
    t12 = evalin('base', 't12');
    t13 = evalin('base', 't13');
    t23 = evalin('base', 't23');
    % 计算Δt12, Δt13, Δt23
    delta_t12 = delta_t(t12,tau_ij_obs(1));
    delta_t13 = delta_t(t13,tau_ij_obs(2));
    delta_t23 = delta_t(t23,tau_ij_obs(3));

    % 计算目标函数，即式(4)
    F = (delta_t12^2 + delta_t13^2 + delta_t23^2) / 75;
end


function shifted_signal = shift_signal(signal, shift_amount)

    % 使用 circshift 进行平移
    shifted_signal = circshift(signal, shift_amount);
    % 如果是向左平移，右侧补零；如果是向右平移，左侧补零
    if shift_amount < 0
        shifted_signal(end+shift_amount+1:end) = 0;
    else
        shifted_signal(1:shift_amount) = 0;
    end
    
end


function tau = showfitted(data)
    % 从数据中找到y的最大值及其索引
    [~, max_index] = max(data(:, 2));
    % 获取最大值周围的3个点的索引
    if max_index > length(data(:,2))-3 || max_index <3
        tau = 20/0.299552816 + 1;
    else
        fit_range = [ -3,-2,-1, 0, 1, 2, 3] + max_index;
        % 获取10个点的索引和对应的值
        fit_indices = data(fit_range, 1);
        fit_values = data(fit_range, 2);
        % 进行抛物线拟合
        coefficients = polyfit(fit_indices, fit_values, 2);
        % 根据拟合结果计算拟合曲线上的点
        fit_indices_curve = linspace(min(fit_indices), max(fit_indices), 1000);
        fit_values_curve = polyval(coefficients, fit_indices_curve);
        % 绘制原始数据和拟合曲线
%         figure;
%         plot(data(:, 1), data(:, 2))
%         plot(data(:, 1), data(:, 2), 'b', fit_indices_curve, fit_values_curve, 'r--');
%         legend('原始数据', '拟合曲线');
%         xlabel('y的索引');
%         ylabel('y的值');
        [~, max_index_fit] = max(fit_values_curve);
        tau = fit_indices_curve(1,max_index_fit);
    end
    
end


%函数：对主窗口进行上采样
function new_signal = upsampling(original_signal,upsampling_factor)

    % 原信号
    original_x = (1:numel(original_signal))';
    original_y = original_signal;
    % 上采样后的采样点数
    upsampled_length = length(original_x) * upsampling_factor;
    % 上采样后的采样点的 x 坐标
    upsampled_x = linspace(1, length(original_x), upsampled_length);
    % 使用多项式插值对原信号进行上采样
    interpolated_signal = interp1(original_x, original_y, upsampled_x, 'spline');
    new_signal = [upsampled_x; interpolated_signal];
end


%对互相关函数进行上采样
function upsampling_gcc = upsampling_gc(r,lag,upsampling_factor)

    % 上采样后的采样点数
    upsampled_length = length(lag) * upsampling_factor;
    % 上采样后的采样点的 x 坐标
    upsampled_x = linspace(-numel(r)/2, numel(r)/2, upsampled_length);
    % 使用多项式插值对原信号进行上采样
    interpolated_signal = interp1(lag, r, upsampled_x, 'spline');
    upsampling_gcc = [upsampled_x; interpolated_signal]';

end


function windowed_signal = windowsignal(signal)
%     r_length = length(signal);
%    % 使用汉明窗
%    window = hamming(r_length);
%    % 对滤波后的信号应用窗函数
%    windowed_signal = signal .* window; % 信号与窗函数相乘
% 
    X = fft(signal);      %变换到频域加窗
    r_length = length(X);
    window = hamming(r_length);
%     得到的是频域信号
    X_windowed = X .* window;

% %     % 进行逆傅里叶变换得到时域信号
      windowed_signal = ifft(X_windowed);

end




function delay = cal_delay(R_xy)
    r_xy = ifft(R_xy);

% 找到主峰值位置
[~, max_idx] = max(r_xy);

% 计算估计的时间延迟
delay = max_idx / 200e3;

end