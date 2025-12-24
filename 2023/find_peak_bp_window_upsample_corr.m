clear;
% === 全局参数 ===
N = 3;
c = 0.299792458;
fs = 200e6;
fp_start = 20e6; % 通带起始
fp_end = 80e6;   % 通带结束
upsampling_factor = 50;
min_peak_distance = 128;
processing_window_len = 512;
window_type = 'hann'; % 重命名变量以免与 loop 中的 window 混淆

% === 几何与文件参数 ===
signal_length = 3e8; % 总处理长度
r_loction = 3e8;     % 起始绝对位置
d12 = 24.9586;
d13 = 34.9335;
d23 = 24.9675;
angle12 = -110.8477;
angle13 = -65.2405;
angle23 = -19.6541;

% === 分段处理参数 ===
block_size = 2e7;    % 每次处理 2000万点 
overlap = 5000;      % 重叠长度 
step_size = block_size - overlap; % 每次前进的步长

% === 结果文件初始化 ===
file_name = 'results\20230718175104_result_yld_3e8_6e8_window_512_128_阈值4倍标准差_去零飘_'+string(fp_start/1e6)+'_'+string(fp_end/1e6)+'_'+ window_type +'.txt';
fileID = fopen(file_name, 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','peak','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');

% === 阈值计算 ===
noise = read_signal('20230718175104.9180CH3.dat', 1e5, 1e8);
filtered_noise = filter_bp(noise, fp_start, fp_end, 5);
min_peak_height = mean(filtered_noise) + 3 * std(filtered_noise);
clear noise filtered_noise; % 释放内存

% === 准备窗函数 ===
if strcmp(window_type,'hann')
    win = hann(processing_window_len);
elseif strcmp(window_type,'exp_hann')
    win = exp_hanning(processing_window_len);
elseif strcmp(window_type,'blackman')
    win = blackman(processing_window_len);
elseif strcmp(window_type,'gausswin')
    win = gausswin(processing_window_len, 8);
else
    win = [];
end
if ~isempty(win), win = win(:); end

% === 进度条 ===
h = waitbar(0, '开始分段处理...');

% === 分段循环 ===
current_read_offset = 0; % 相对 r_loction 的偏移量

while current_read_offset < signal_length
    
    % 1. 计算当前读取长度
    current_len = min(block_size, signal_length - current_read_offset);
    
    % 如果剩余数据不足以构成一个处理窗口，则停止
    if current_len < processing_window_len
        break;
    end
    
    % 更新进度条
    waitbar(current_read_offset / signal_length, h, ...
        sprintf('正在处理分段: %.1f%% (Offset: %d)', (current_read_offset/signal_length)*100, current_read_offset));
    
    % 2. 读取当前分段数据
    % 注意：read_signal 的第三个参数是绝对位置
    current_abs_loc = r_loction + current_read_offset;
    
    ch1 = read_signal('20230718175104.9180CH1.dat', current_len, current_abs_loc);
    ch2 = read_signal('20230718175104.9180CH2.dat', current_len, current_abs_loc);
    ch3 = read_signal('20230718175104.9180CH3.dat', current_len, current_abs_loc);
    
    % 3. 滤波 (对分段数据滤波)
    processed_ch1_yld = filter_bp(detrend(ch1), fp_start, fp_end, 5);
    processed_ch2_yld = filter_bp(detrend(ch2), fp_start, fp_end, 5);
    processed_ch3_yld = filter_bp(detrend(ch3), fp_start, fp_end, 5);
    
    % 释放原始数据内存
    clear ch1 ch2 ch3;
    
    % 4. 寻找峰值 (在当前分段内)
    [pks, locs] = findpeaks(processed_ch1_yld, ...
        'MinPeakHeight', min_peak_height, ...
        'MinPeakDistance', min_peak_distance);
    
    num_peaks = length(locs);
    
    % 5. 遍历处理峰值
    for i = 1:num_peaks
        center_loc = locs(i); % 分段内的局部坐标
        
        % --- 关键：分段边界与重叠处理逻辑 ---
        
        % A. 截取范围计算
        start_idx = center_loc - floor(processing_window_len / 2);
        end_idx = center_loc + ceil(processing_window_len / 2) - 1;
        
        % B. 物理边界检查 (分段开头)
        if start_idx < 1
            continue; 
        end
        
        % C. 物理边界检查 (分段结尾)
        if end_idx > current_len
            continue;
        end
        
        % D. 重叠区域去重
        % 如果这不是最后一段数据，且峰值位于当前段的末尾重叠区 (Overlap Region)，
        % 则跳过该峰值。这个峰值会在下一段数据的开头（非边缘区域）被再次检测到并处理。
        % 这样做是为了保证边缘处的峰值能获得完整的上下文，且避免滤波边缘效应。
        is_last_block = (current_read_offset + current_len >= signal_length);
        if ~is_last_block && (center_loc > (current_len - overlap))
            continue; 
        end
        
        % --- 信号截取与处理 ---
        segment_ch1 = processed_ch1_yld(start_idx:end_idx);
        segment_ch2 = processed_ch2_yld(start_idx:end_idx);
        segment_ch3 = processed_ch3_yld(start_idx:end_idx);
        
        % 加窗
        if isempty(win)
            % 假设 template_norm 在代码其他地方定义了，或者不需要
             % windowed_segment_ch1 = conv_window(segment_ch1, template_norm);
             % 如果没有 win 且没有 template_norm 逻辑，保持原样：
             windowed_segment_ch1 = segment_ch1; 
             windowed_segment_ch2 = segment_ch2;
             windowed_segment_ch3 = segment_ch3;
        else
            windowed_segment_ch1 = segment_ch1 .* win;
            windowed_segment_ch2 = segment_ch2 .* win;
            windowed_segment_ch3 = segment_ch3 .* win;
        end
        
        % 上采样
        [ch1_up, ch2_up, ch3_up] = deal(...
            upsampling(windowed_segment_ch1, upsampling_factor)', ...
            upsampling(windowed_segment_ch2, upsampling_factor)', ...
            upsampling(windowed_segment_ch3, upsampling_factor)');
        
        ch1_upsp = ch1_up(:,2);
        ch2_upsp = ch2_up(:,2);
        ch3_upsp = ch3_up(:,2);
        
        % 互相关
        [r12_gcc, lags12_gcc] = xcorr(ch1_upsp, ch2_upsp, 'normalized');
        [r13_gcc, lags13_gcc] = xcorr(ch1_upsp, ch3_upsp, 'normalized');
        [r23_gcc, lags23_gcc] = xcorr(ch2_upsp, ch3_upsp, 'normalized');
        
        R12_gcc = max(r12_gcc);
        R13_gcc = max(r13_gcc);
        R23_gcc = max(r23_gcc);
        
        t12_gcc = cal_tau(r12_gcc, lags12_gcc');
        t13_gcc = cal_tau(r13_gcc, lags13_gcc');
        t23_gcc = cal_tau(r23_gcc, lags23_gcc');
        
        % 还原时延
        t12 = t12_gcc * 5 / upsampling_factor;
        t13 = t13_gcc * 5 / upsampling_factor;
        t23 = t23_gcc * 5 / upsampling_factor;
        
        % 几何解算初值
        cos_beta_0 = ((c*t13*d12*sind(angle12))-(c*t12*sind(angle13)*d13))/(d13*d12*sind(angle12-angle13));
        cos_alpha_0 = ((c*t12)/d12-cos_beta_0*cosd(angle12))/sind(angle12);
        
        if abs(cos_beta_0) > 1 || abs(cos_alpha_0) > 1
            continue;
        end
        
        x0 = [cos_alpha_0, cos_beta_0];
        
        % 优化求解
        % 关闭 Display 以提高大量循环的速度
        options = optimoptions('lsqnonlin', 'MaxIter', 1000, 'TolFun', 1e-6, 'Display', 'off');
        x = lsqnonlin(@(x) objective(x, t12, t13, t23,'yld'), x0, [-1 -1], [1 1], options);
        
        cos_alpha_opt = x(1);
        cos_beta_opt = x(2);
        
        if abs(cos_alpha_opt) > 1 || abs(cos_beta_opt) > 1
            continue;
        end
        
        Az = atan2(cos_alpha_opt, cos_beta_opt);
        if abs(cos_beta_opt/cos(Az)) > 1
            continue;
        end
        
        El = acos(cos_beta_opt/cos(Az));
        Az_deg = rad2deg(Az);
        El_deg = rad2deg(El);
        if Az_deg < 0
            Az_deg = Az_deg + 360;
        end
        
        t123 = t12 + t23 - t13;
        Rcorr = (R12_gcc + R13_gcc + R23_gcc) / 3;
        
        % 写入结果
        % 全局坐标 = 基准 + 当前段偏移 + 窗口起始局部坐标
        global_start_loc = r_loction + current_read_offset + start_idx;
        
        fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
            global_start_loc, processing_window_len, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr, t123);
            
    end % end peak loop
    
    % 更新偏移量，准备读取下一段
    % 步长 = 块大小 - 重叠量
    current_read_offset = current_read_offset + step_size;
    
end % end while loop

fclose(fileID);
close(h);
fprintf('处理完成。\n');

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


function tau = cal_tau(R, lag)
% 从数据中找到y的最大值及其索引
[~, max_index] = max(R);
tau = lag(max_index,1);
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



function w = exp_hanning(n, alpha)
% 指数加权汉宁窗：通过指数因子强化边缘衰减
% 输入：n - 窗长；alpha - 陡峭度参数（>0，值越大边缘越陡）
% 输出：w - 指数加权汉宁窗（归一化至最大值为1）

if nargin < 2
    alpha = 3;  % 默认陡峭度参数
end

% 生成0到1的归一化索引
k = 0:n-1;
% 标准汉宁窗
hann_win = 0.5 - 0.5 * cos(2*pi*k/(n-1));
% 指数因子：中心权重为1，向边缘快速衰减
exp_factor = exp(-alpha * (abs(k - (n-1)/2) ./ ((n-1)/2)).^2);
% 组合并归一化
w = hann_win .* exp_factor;
w = w / max(w);
end


%% 设计巴特沃斯带通滤波器
function filtered_signal = filter_bp(signal,f1,f2,order)
    Fs = 200e6;
    fn = Fs/2;
    Wn = [f1 f2]/fn;
    [b,a] = butter(order,Wn); 
    filtered_signal = filtfilt(b,a,signal);

end