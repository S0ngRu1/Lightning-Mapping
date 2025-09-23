%% ========================================================================
%  二维定位算法
% =========================================================================
clear; clc; close all;
fp_start = 30e6; % 通带起始
fp_end = 80e6;   % 通带结束
N = 3; c = 0.299792458; fs = 200e6; step = 1e4; upsampling_factor = 1;
start_signal_loc = 3.65e8; end_signal_loc = 3.7e8;
all_start_signal_loc = start_signal_loc:step:end_signal_loc;
%引雷点
d12 = 24.9586; d13 = 34.9335; d23 = 24.9675;
angle12 = -110.8477; angle13 = -65.2405; angle23 = -19.6541;
%引雷点阈值
noise = read_signal('20240613171023.1120CH1.dat',1e5,1e8);
filtered_noise = filter_bp(noise,fp_start,fp_end,5);
noise_std = std(filtered_noise);
threshold_factor = 3;      % find_pulses_advanced 的阈值因子
merge_gap_samples = 10;   % 脉冲融合的间隙阈值
pulses_per_group = 4; % 定义每组包含n个脉冲
% --- 文件写入准备 ---
file_name = '20240613171023_result_yld_3.65e8_3.7e8_find_pulse_阈值3倍标准差_去零飘_'+string(fp_start/1e6)+'_'+string(fp_end/1e6)+'.txt';
fileID = fopen(file_name, 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','Pulse_Len','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');
num_total_blocks = numel(all_start_signal_loc)-1;

h_overall = waitbar(0, '正在初始化处理...', 'Name', '整体处理进度');
for j = 1:num_total_blocks

    progress = j / num_total_blocks;
    waitbar(progress, h_overall, sprintf('整体进度: 正在处理信号块 %d / %d', j, num_total_blocks));
    current_block_start = all_start_signal_loc(j);
    % --- 1. 读取当前处理块的完整信号 ---
    %     引雷点
    ch1 = read_signal('20240613171023.1120CH1.dat', step, current_block_start);
    ch2 = read_signal('20240613171023.1120CH2.dat', step, current_block_start);
    ch3 = read_signal('20240613171023.1120CH3.dat', step, current_block_start);
    filtered_signal1 = filter_bp(detrend(ch1),fp_start,fp_end,5);
    filtered_signal2 = filter_bp(detrend(ch2),fp_start,fp_end,5);
    filtered_signal3 = filter_bp(detrend(ch3),fp_start,fp_end,5);
    pulse_catalog_in_block = find_pulses_advanced(filtered_signal1, noise_std, fs, threshold_factor, merge_gap_samples);

    if isempty(pulse_catalog_in_block)
        fprintf('      在本块内未找到有效脉冲，跳过。\n');
        continue;
    end

    num_pulses_in_block = numel(pulse_catalog_in_block);
    fprintf('      在本块内找到 %d 个精确脉冲，开始处理...\n', num_pulses_in_block);

    % --- 遍历当前块内找到的【脉冲目录】 ---
    for pi = 1:pulses_per_group: num_pulses_in_block
        start_group_idx = pi;
        end_group_idx = min(pi + pulses_per_group - 1, num_pulses_in_block); % 防止末尾不足一组时越界
        pulse_group = pulse_catalog_in_block(start_group_idx:end_group_idx);
        % 如果组内脉冲太少，可以跳过
        if numel(pulse_group) < 2
            continue;
        end


        % **  确定该组的总起点和总终点 **
        group_start_idx = pulse_group(1).start_idx;
        group_end_idx   = pulse_group(end).end_idx;

        % --- 3. 根据总起止点，截取包含整组脉冲的长信号片段 ---
        signal1 = filtered_signal1(group_start_idx : group_end_idx);
        signal2 = filtered_signal2(group_start_idx : group_end_idx);
        signal3 = filtered_signal3(group_start_idx : group_end_idx);

        pulse_len = length(signal1); % 这是整个脉冲组的总长度

        % ** 长度对齐 (零填充)，以保证互相关等操作正常进行 **
        max_len = max([length(signal1), length(signal2), length(signal3)]);
        signal1_padded = [signal1; zeros(max_len - length(signal1), 1)];
        signal2_padded = [signal2; zeros(max_len - length(signal2), 1)];
        signal3_padded = [signal3; zeros(max_len - length(signal3), 1)];

        % --- 后续的所有处理都使用填充对齐后的信号 ---
        [ch1_new, ch2_new, ch3_new] = deal(...
            real(windowsignal(signal1_padded)), ...
            real(windowsignal(signal2_padded)), ...
            real(windowsignal(signal3_padded)));


        % 上采样
        [ch1_up, ch2_up, ch3_up] = deal(...
            upsampling(ch1_new, upsampling_factor)', ...
            upsampling(ch2_new, upsampling_factor)', ...
            upsampling(ch3_new, upsampling_factor)');
        ch1_upsp = ch1_up(:,2);
        ch2_upsp = ch2_up(:,2);
        ch3_upsp = ch3_up(:,2);

        %互相关
        [r12_gcc,lags12_gcc] = xcorr(ch1_upsp ,ch2_upsp,'normalized');
        [r13_gcc,lags13_gcc] = xcorr(ch1_upsp ,ch3_upsp,'normalized');
        [r23_gcc,lags23_gcc] = xcorr(ch2_upsp ,ch3_upsp,'normalized');
        R12_gcc = max(r12_gcc);
        R13_gcc = max(r13_gcc);
        R23_gcc = max(r23_gcc);
        t12_gcc = cal_tau(r12_gcc,lags12_gcc');
        t13_gcc = cal_tau(r13_gcc,lags13_gcc');
        t23_gcc = cal_tau(r23_gcc,lags23_gcc');

        %
        %         %从化局
        %         t12 = t12_gcc *5/upsampling_factor;
        %         t13 = t13_gcc *5/upsampling_factor+1.600061;
        %         t23 = t23_gcc *5/upsampling_factor+1.600061;


        %引雷场
        t12 = t12_gcc *5/upsampling_factor;
        t13 = t13_gcc *5/upsampling_factor;
        t23 = t23_gcc *5/upsampling_factor;

        cos_beta_0 =((c*t13*d12*sind(angle12))-(c*t12*sind(angle13)*d13))/(d13*d12*sind(angle12-angle13)) ;
        cos_alpha_0 = ((c*t12)/d12-cos_beta_0*cosd(angle12))/sind(angle12);
        if abs(cos_beta_0)>1 || abs(cos_alpha_0)>1
            continue;
        end
        x0 = [cos_alpha_0,cos_beta_0];
        % 调用lsqnonlin函数进行优化
        options = optimoptions('lsqnonlin', 'MaxIter', 1000, 'TolFun', 1e-6);
        x = lsqnonlin(@(x) objective(x, t12, t13, t23,'yld'), x0, [-1 -1],[1 1], options);
        % 输出最优的cos(α)和cos(β)值
        cos_alpha_opt = x(1);
        cos_beta_opt = x(2);
        if abs(cos_alpha_opt)>1 || abs(cos_beta_opt)>1
            continue;
        end
        Az = atan2( cos_alpha_opt,cos_beta_opt);
        if abs(cos_beta_opt/cos(Az)) > 1
            continue;
        end
        El = acos( cos_beta_opt/cos(Az) );
        % 将弧度转换为角度
        Az_deg = rad2deg(Az);
        El_deg = rad2deg(El);
        if Az_deg < 0
            Az_deg = Az_deg + 360;
        end

        t123 = t12 + t23 - t13;
        Rcorr = (R12_gcc + R13_gcc + R23_gcc) / 3;
        % --- 写入结果 ---
        ts_ns = 1/fs*1e9;
        absolute_sample_location = current_block_start + group_start_idx - 1;
        fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
            absolute_sample_location, pulse_len, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
    end

end
close(h_overall);
fclose(fileID);


result1 = readtable(file_name);


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


function windowed_signal = windowsignal(signal)
%     r_length = length(signal);
%    % 使用汉明窗
%    window = hamming(r_length);
%    % 对滤波后的信号应用窗函数
%    windowed_signal = signal .* window; % 信号与窗函数相乘

    X = fft(signal);      %变换到频域加窗
    r_length = length(X);
    window = hamming(r_length);
%     得到的是频域信号
    X_windowed = X .* window;

% %     % 进行逆傅里叶变换得到时域信号
      windowed_signal = ifft(X_windowed);

end

