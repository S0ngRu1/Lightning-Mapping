% === 全局参数设置 ===
total_signal_length = 2e8; % 总信号长度
r_loction_base = 3e8;      % 初始的绝对位置基准

% --- 分段参数 ---
block_size = 1e7;          % 每次处理 1千万点 (约占用几百MB内存)
overlap = 5000;            % 重叠区域 (必须 > bigwindows_length)
step_size = block_size - overlap; % 每次循环前进的步长

% --- 算法参数 ---
N = 3;
d12 = 24.96;
d13 = 34.93;
d23 = 24.98;
c = 0.299552816;
fs = 200e6;
upsampling_factor = 50;
window_length = 1024;
bigwindows_length = window_length + 100;
window = window_length * upsampling_factor;
msw_length = 50;
angle12 = -110.85;
angle13 = -65.24;
angle23 = -19.65;

% 打开结果文件 (只打开一次)
result_filename = '20230718175104_result_3e8_5e8_window_1024_256_corr_msw_mean.txt';
fileID = fopen(result_filename, 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','peak','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');

% === 分段处理循环 ===
current_start_idx = 0; % 当前读取的起始点 (0-based)

while current_start_idx < total_signal_length
    
    % 1. 计算当前段的读取长度
    % 如果是最后一段，读取长度可能会小于 block_size
    current_read_len = min(block_size, total_signal_length - current_start_idx);
    
    if current_read_len < bigwindows_length
        break; % 如果剩余数据不够一个窗口，停止
    end
    
    fprintf('正在处理分段: 起始点 %.0f, 长度 %.0f\n', current_start_idx, current_read_len);

    % 2. 读取信号 (注意：这里需要你的 read_signal 支持偏移量)
    % 假设 read_signal 现在的用法是: read_signal(文件名, 读取长度, 起始偏移量)
    % 如果你的 read_signal 不支持第三个参数，你需要修改它使用 fseek 定位
    
    % 这里的 r_loction 参数仅仅传给读取函数用于定位文件指针
    % 我们假设数据文件按顺序存储
    seek_offset = current_start_idx; 
    
    % [重要] 请确保你的 read_signal 函数内部使用了 fseek(fid, seek_offset * bytes_per_sample, 'bof')
    ch1 = read_signal('20230718175104.9180CH1.dat', current_read_len, seek_offset);
    ch2 = read_signal('20230718175104.9180CH2.dat', current_read_len, seek_offset);
    ch3 = read_signal('20230718175104.9180CH3.dat', current_read_len, seek_offset);

    % 3. 滤波 (对当前分段进行滤波)
    % 注意：分段滤波在边缘会有瞬态效应，但由于我们有 overlap，下次循环会覆盖边缘
    filtered_signal1 = filter_bp(ch1, 20e6 ,80e6 ,5);
    filtered_signal2 = filter_bp(ch2, 20e6 ,80e6 ,5);
    filtered_signal3 = filter_bp(ch3, 20e6 ,80e6 ,5);

    % 4. 计算统计量 (基于当前分段)
    mean_abs_signal1 = mean(abs(filtered_signal1));
    % std_signal1 = std(filtered_signal1); % 未使用，注释掉节省计算
    threshold = 0.05 * mean_abs_signal1;

    % 5. 平滑
    smoothed_signal1 = movmean(filtered_signal1, 10);
    smoothed_signal2 = movmean(filtered_signal2, 10);
    smoothed_signal3 = movmean(filtered_signal3, 10);

    % 6. 寻找峰值 (局部坐标)
    [peaks, locs] = findpeaks(smoothed_signal1, 'MinPeakHeight', threshold, 'MinPeakDistance', window_length/4);

    num_peaks = numel(peaks);
    for pi = 1:num_peaks
        idx = locs(pi); % 这是一个局部索引 (1 到 current_read_len)
        
        % --- 边界检查与重叠去重 ---
        % 1. 窗口越界检查：确保峰值周围有足够的数据截取窗口
        if idx - (bigwindows_length / 2 - 1) <= 0 || idx + (bigwindows_length / 2) > length(smoothed_signal1)
            continue;
        end
        
        % 2. [关键] 重叠区域去重处理
        % 我们只处理当前段 "独有" 的部分，或者是第一段。
        % 如果峰值落在 overlap 区域（即当前段的末尾），我们跳过它，
        % 因为它会在下一段的开头作为非边缘信号被再次处理。
        % 这样可以避免同一个峰值被处理两次，也可以避免处理位于段边缘滤波效果不好的峰值。
        
        % 有效区结束点：如果是最后一段，则处理到最后；否则保留 overlap 区域给下一段
        valid_end_idx = current_read_len; 
        if (current_start_idx + current_read_len) < total_signal_length
            valid_end_idx = current_read_len - overlap; 
        end
        
        % 如果是第一段之后的数据，为了防止重复，也要忽略开头的一小段(由上一段覆盖)
        % 但通常 filter_bp 在开头会有暂态，保留 overlap 给下一段开头是比较好的策略。
        % 简单策略：只处理 idx <= valid_end_idx
        if idx > valid_end_idx
            continue; 
        end
        
        % ---------------------------------------------------------
        % 下面是原本的核心处理逻辑 (几乎不需要改动，除了全局坐标计算)
        % ---------------------------------------------------------
        
        % 截取窗口信号
        win_start = idx-(bigwindows_length/2-1);
        win_end = idx+(bigwindows_length/2);
        
        [signal1, signal2, signal3] = deal(...
            smoothed_signal1(win_start:win_end), ...
            smoothed_signal2(win_start:win_end), ...
            smoothed_signal3(win_start:win_end));
            
        % ... (此处省略中间未变动的去直流、上采样、互相关逻辑，保持原样) ...
        [ch1_new, ch2_new, ch3_new] = deal(real(windowsignal(detrend(signal1))), real(windowsignal(detrend(signal2))), real(windowsignal(detrend(signal3))));
        [ch1_up, ch2_up, ch3_up] = deal(upsampling(ch1_new, upsampling_factor)', upsampling(ch2_new, upsampling_factor)', upsampling(ch3_new, upsampling_factor)');
        ch1_upsp = ch1_up(:,2); ch2_upsp = ch2_up(:,2); ch3_upsp = ch3_up(:,2);
        
        % 截取中心部分做互相关
        center_idx_up = bigwindows_length*upsampling_factor/2;
        win_half_up = window/2;
        slice_idx = center_idx_up-(win_half_up-1) : center_idx_up+win_half_up;
        
        ch1_new_corr = ch1_upsp(slice_idx);
        ch2_new_corr = ch2_upsp(slice_idx);
        ch3_new_corr = ch3_upsp(slice_idx);
        
        [r12_gcc,lags12_gcc] = xcorr(ch1_new_corr,ch2_new_corr,'normalized');
        [r13_gcc,lags13_gcc] = xcorr(ch1_new_corr,ch3_new_corr,'normalized');
        [r23_gcc,lags23_gcc] = xcorr(ch2_new_corr,ch3_new_corr,'normalized');
        
        R12_gcc = max(r12_gcc); R13_gcc = max(r13_gcc); R23_gcc = max(r23_gcc);
        t12_gcc = cal_tau(r12_gcc,lags12_gcc');
        t13_gcc = cal_tau(r13_gcc,lags13_gcc');
        t23_gcc = cal_tau(r23_gcc,lags23_gcc');
        
        ismsw = 0;
        
        % 逻辑分支：只有相关性高才计算 DOA
        if R12_gcc > 0.8 && R13_gcc > 0.8
            shifted_ch1 = ch1_new_corr;
            shifted_ch2 = shift_signal(ch2_upsp, t12_gcc); shifted_ch2 = shifted_ch2(slice_idx);
            shifted_ch3 = shift_signal(ch3_upsp, t13_gcc); shifted_ch3 = shifted_ch3(slice_idx);
            
            peaks1 = find_peaks(shifted_ch1, mean(shifted_ch1));
            peaks2 = find_peaks(shifted_ch2, mean(shifted_ch2));
            peaks3 = find_peaks(shifted_ch3, mean(shifted_ch3));
            matched_peaks_x = match_peaks(peaks1, peaks2, peaks3);
            
            t12s = []; t13s = []; t23s = [];
            
            if ~isempty(matched_peaks_x)
                for i = 1 : size(matched_peaks_x, 1)
                     %微尺度 (注意：这里需确保 msw_signal 函数内部实现没问题)
                    ch1_msw = msw_signal(shifted_ch1 , matched_peaks_x(i,1) ,msw_length,window);
                    ch2_msw = msw_signal(shifted_ch2 , matched_peaks_x(i,2) ,msw_length,window);
                    ch3_msw = msw_signal(shifted_ch3 , matched_peaks_x(i,3) ,msw_length,window);
                    
                    if numel(ch1_msw)~= msw_length*2 || numel(ch2_msw)~= msw_length*2 || numel(ch3_msw)~= msw_length*2
                        continue;
                    end
                    
                    [R12_msw,lags12_msw] = xcorr(ch1_msw,ch2_msw,'normalized');
                    [R13_msw,lags13_msw] = xcorr(ch1_msw,ch3_msw,'normalized');
                    [R23_msw,lags23_msw] = xcorr(ch2_msw,ch3_msw,'normalized');
                    
                    if max(R12_msw) > 0.8 && max(R13_msw) > 0.8 && max(R23_msw) > 0.8
                         t12s = [t12s cal_tau(R12_msw,lags12_msw')];
                         t13s = [t13s cal_tau(R13_msw,lags13_msw')];
                         t23s = [t23s cal_tau(R23_msw,lags23_msw')];
                    end
                end
                
                if ~isempty(t12s) && ~isempty(t13s) && ~isempty(t23s)
                    t12 = (t12_gcc + mean(t12s))*5/upsampling_factor;
                    t13 = (t13_gcc + mean(t13s))*5/upsampling_factor;
                    t23 = (t23_gcc + mean(t23s))*5/upsampling_factor;
                    
                    % --- DOA 计算 (封装成函数或保持原样) ---
                    [Az_deg, El_deg, cos_alpha_opt, cos_beta_opt, valid] = calculate_doa(t12, t13, d12, d13, angle12, angle13, c);
                    
                    if valid
                        t123 = t12 + t23 - t13;
                        Rcorr = (R12_gcc + R13_gcc + R23_gcc) / 3;
                        
                        % [关键] 全局坐标计算
                        % 全局位置 = 基准 + 当前段起始偏移 + 段内索引
                        global_loc = r_loction_base + current_start_idx + idx - window/upsampling_factor/2; 
                        
                        fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
                            global_loc, window/upsampling_factor/2, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr, t123);
                        ismsw = ismsw + 1;
                    end
                end
            end
        end
        
        % 如果 MSW 没有结果，退回 GCC 结果
        if ismsw == 0
            t12 = t12_gcc *5/upsampling_factor;
            t13 = t13_gcc *5/upsampling_factor;
            t23 = t23_gcc *5/upsampling_factor;
            
            [Az_deg, El_deg, cos_alpha_opt, cos_beta_opt, valid] = calculate_doa(t12, t13, d12, d13, angle12, angle13, c);
            
            if valid
                t123 = t12 + t23 - t13;
                Rcorr = (R12_gcc + R13_gcc + R23_gcc) / 3;
                
                % [关键] 全局坐标计算
                global_loc = r_loction_base + current_start_idx + idx - window/100;
                
                fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
                    global_loc, window/100, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr, t123);
            end
        end
        
    end % end of peaks loop
    
    % 更新起始点，准备读取下一段
    current_start_idx = current_start_idx + step_size;
    
end % end of while loop

fclose(fileID);
fprintf('处理完成。\n');


function delta_t = delta_t(tij,tij_obs)
    delta_t = tij - tij_obs;
end




function tau = cal_tau(R, lag)
    % 从数据中找到y的最大值及其索引
    [~, max_index] = max(R);
    tau = lag(max_index,1);
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


function mswed_signal = msw_signal(signal , peak_x ,length,win)
      % 找到峰值的 x 值在信号中的索引
    left_idx = max(peak_x - length+1, 1);  % 确定左边界的索引
    right_idx = min(peak_x + length, win);  % 确定右边界的索引
    mswed_signal = signal(left_idx:right_idx);  % 提取以中心 x 值为中心的左右40个采样点

end


% 定义目标函数 (正确版本)
function F = objective(x, t12_meas, t13_meas, t23_meas, type)
% 提取待优化的变量
cos_alpha = x(1);
cos_beta = x(2);

% 计算τij的理论值 τ_model (我将 obs 改为 model，语义更清晰)
tau_model = calculate_tau_obs(cos_alpha, cos_beta, type);

% t12, t13, t23 是测量的时延 (measurement)
% tau_model(1), tau_model(2), tau_model(3) 是根据当前 x 计算出的理论时延

% 计算残差向量
residual12 = t12_meas - tau_model(1);
residual13 = t13_meas - tau_model(2);
residual23 = t23_meas - tau_model(3);

% 返回残差向量 F
% lsqnonlin 会自动最小化 sum(F.^2)
F = [residual12; residual13; residual23];
end


% 定义计算τij的理想值τ_ij^obs的函数
function tau_ij_obs = calculate_tau_obs(cos_alpha, cos_beta, type)
% 初始化输出变量
tau_ij_obs = zeros(1, 3);

% 根据 type 参数选择不同的参数集
if strcmp(type, 'chj') % 从化局
    angle12 = -2.8381;
    angle13 = 50.3964;
    angle23 = 120.6568;
    d12 = 41.6496;
    d13 = 36.9015;
    d23 = 35.4481;
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



function signal = read_signal(signal_path, r_length,r_loction)
fid  = fopen(signal_path,'r');%读取数据的位置

%使用fseek函数将文件指针移动到指定位置，以便读取数据。
%这里指定移动位置为r_location，表示移动到指定位置开始读取数据。
fseek(fid,r_loction*2,'bof');
%使用fread函数从文件中读取数据，读取的数据长度为r_length，数据以int16格式读取。
%将读取到的数据分别保存到变量ch_1、ch_2和ch_3中。
signal = fread(fid,r_length,'int16');
%关闭所有文件
fclose(fid);
end
