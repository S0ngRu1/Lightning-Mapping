clear;
N = 3;
c = 0.299792458;
fs = 200e6;
fp_start = 30e6; % 通带起始
fp_end = 80e6;   % 通带结束
upsampling_factor = 50;
window = 'hann';
% window = 'exp_hann';
% window = 'blackman';
% window = 'gausswin';
% 从化局
% angle12 = -2.8381;
% angle13 = 50.3964;
% angle23 = 120.6568;
% d12 = 41.6496;
% d13 = 36.9015;
% d23 = 35.4481;
% signal_length = 2.1e8;
% r_loction = 3.9e8;
% ch1 = read_signal('..\\2024 822 85933.651462CH1.dat',signal_length,r_loction);
% ch2 = read_signal('..\\2024 822 85933.651462CH2.dat',signal_length,r_loction);
% ch3 = read_signal('..\\2024 822 85933.651462CH3.dat',signal_length,r_loction+215/5);
%引雷点
signal_length = 2e7;
r_loction = 3.8e8;
d12 = 24.9586;
d13 = 34.9335;
d23 = 24.9675;
angle12 = -110.8477;
angle13 = -65.2405;
angle23 = -19.6541;
ch1 = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction);
ch2 = read_signal('..\\20240822165932.6610CH2.dat',signal_length,r_loction);
ch3 = read_signal('..\\20240822165932.6610CH3.dat',signal_length,r_loction);

% 设计一个8阶巴特沃斯带通滤波器
% order = 8;
% [b, a] = butter(order, [fp_start, fp_end] / (fs/2), 'bandpass');
% % 对信号进行零相位滤波
% processed_ch1_yld = filtfilt(b, a, ch1);
% processed_ch2_yld = filtfilt(b, a, ch2);
% processed_ch3_yld = filtfilt(b, a, ch3);

processed_ch1_yld = filter_bp(detrend(ch1),fp_start,fp_end,5);
processed_ch2_yld = filter_bp(detrend(ch2),fp_start,fp_end,5);
processed_ch3_yld = filter_bp(detrend(ch3),fp_start,fp_end,5);


file_name = 'result_yld_3.8e8_4e8_window_512_128_阈值60_去零飘1_'+string(fp_start/1e6)+'_'+string(fp_end/1e6)+'_'+ window +'.txt';
% 打开一个文本文件用于写入运行结果
fileID = fopen(file_name, 'w');
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','peak','t12', 't13', 't23', 'cos_alpha_opt', 'cos_beta_opt','Azimuth', 'Elevation', 'Rcorr', 't123');

% %引雷点阈值
% noise = read_signal('..\\20240822165932.6610CH1.dat',1e8,1e8);
% filtered_noise = filter_bp(noise,fp_start,fp_end,5);
% min_peak_height = mean(filtered_noise)+5*std(filtered_noise);
min_peak_height = 60;
min_peak_distance = 128;
% 以峰值为中心，进行处理的信号片段的总长度
processing_window_len =1024;

% 寻找能量峰值的位置
[pks, locs] = findpeaks(processed_ch1_yld, ...
    'MinPeakHeight', min_peak_height, ...
    'MinPeakDistance', min_peak_distance);
if strcmp(window,'hann')
    % 创建一个汉宁窗
    win = hann(processing_window_len);
    win = win(:); % 确保是列向量
elseif strcmp(window,'exp_hann')
    % 创建一个指数汉宁窗
    win = exp_hanning(processing_window_len);
    win = win(:); % 确保是列向量
elseif strcmp(window,'blackman')
    % 创建一个布莱克曼窗
    win = blackman(processing_window_len);
    win = win(:); % 确保是列向量
elseif strcmp(window,'gausswin')
    % 创建一个高斯窗
    win = gausswin(processing_window_len, 8);
    win = win(:); % 确保是列向量
else
    win = [];
end

% 创建进度条
h = waitbar(0, '正在处理峰值...');
num_peaks = length(locs);
for i = 1:num_peaks
    center_loc = locs(i);
    waitbar(i / num_peaks, h, sprintf('正在处理峰值 %d/%d', i, num_peaks));
    % 计算截取片段的起始和结束索引
    start_idx = center_loc - floor(processing_window_len / 2);
    end_idx = center_loc + ceil(processing_window_len / 2) - 1;

    % 边界检查，防止索引越界
    if start_idx < 1 || end_idx > length(processed_ch1_yld)
        fprintf('事件 %d (位置 %d) 过于靠近信号边界，已跳过。\n', i, center_loc);
        continue;
    end

    % 从三个通道截取信号片段
    segment_ch1 = processed_ch1_yld(start_idx:end_idx);
    segment_ch2 = processed_ch2_yld(start_idx:end_idx);
    segment_ch3 = processed_ch3_yld(start_idx:end_idx);
    % 对截取的片段应用窗函数
    if isempty(win)
        windowed_segment_ch1 = conv_window(segment_ch1, template_norm);
        windowed_segment_ch2 = conv_window(segment_ch2, template_norm);
        windowed_segment_ch3 = conv_window(segment_ch3, template_norm);
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


    %从化局
    %         t12 = t12_gcc *0.1;
    %         t13 = t13_gcc *0.1+1.600061;
    %         t23 = t23_gcc *0.1+1.600061;
    %         t13 = t13_gcc *0.1;
    %         t23 = t23_gcc *0.1;

    %     %引雷场
    t12 = t12_gcc *0.1;
    t13 = t13_gcc *0.1;
    t23 = t23_gcc *0.1;

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

    % 写入计算后的数据
    fprintf(fileID, '%-13d%-15d%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
        r_loction+start_idx,processing_window_len, t12, t13, t23, cos_alpha_opt, cos_beta_opt, Az_deg, El_deg, Rcorr,t123);
end
% 关闭文件
fclose(fileID);
% 关闭进度条
close(h);

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


function plot_signal_spectrum(signal)
% Plotting the signal spectrum时域信号的频谱图
fs = 200;
fft_signal = fft(signal);
n = length(fft_signal);
x = (0:n/2-1) * (fs/n);
figure
plot(x, 2.0 / n * abs(fft_signal(1:n/2)))

xlabel('Frequency (MHz)')
ylabel('Amplitude')
grid on
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



%% 设计巴特沃斯带通滤波器
function filtered_signal = filter_bp(signal,f1,f2,order)
    Fs = 200e6;
    fn = Fs/2;
    Wn = [f1 f2]/fn;
    [b,a] = butter(order,Wn); 
    filtered_signal = filtfilt(b,a,signal);

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
    