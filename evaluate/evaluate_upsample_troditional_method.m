clear; clc; close all;

%% === 1. 参数设置 ===
% --- 文件路径 ---
% 1. 第一列：低倍上采样 (Upsample = 2)
file_col1 = '..\2024\2d\results\20240822165932_loop_result_yld_3.92e8_9e6_window_1024_256_ErrCalc_upsample_2.txt';

% 2. 第二列：高倍上采样 (Upsample = 50)
file_col2 = '..\2024\2d\results\20240822165932_loop_result_yld_3.92e8_9e6_window_1024_256_ErrCalc_upsample_50.txt';

% 原始波形文件
raw_file_ch1 = '..\2024\20240822165932.6610CH1.dat'; 
has_waveform = isfile(raw_file_ch1); 

% --- 通用参数 ---
fs = 200e6;              % 采样率
Start_loc_Base = 3.92e8; % 基准采样点位置 (注意：根据您的文件名，这里似乎应该是 3.92e8)
T_us_max = 45000;        % 最大显示时间 (us) - 根据文件长度 9e6/200 ~ 45000us
Downsample_Factor = 100; % 波形降采样倍数

%% === 2. 数据读取与预处理 ===
% 统一筛选参数
t123_th = 1.0;
rcorr_th = 0.5;

% read_and_filter(文件名, 基准点, 闭合差阈值, Rcorr阈值, 最大时间)

% Col 1: Low Upsampling
data_1 = read_and_filter(file_col1, Start_loc_Base, t123_th, rcorr_th, T_us_max);

% Col 2: High Upsampling
data_2 = read_and_filter(file_col2, Start_loc_Base, t123_th, rcorr_th, T_us_max);

% --- 读取原始波形 ---
if has_waveform
    read_len = ceil(T_us_max * 1e-6 * fs); 
    raw_sig = read_signal(raw_file_ch1, read_len, Start_loc_Base);
    t_full = (0:length(raw_sig)-1) / fs * 1e6; 
    
    % 降采样
    raw_sig = raw_sig(1:Downsample_Factor:end);
    t_raw = t_full(1:Downsample_Factor:end); 
    raw_sig = filter_bp(raw_sig,30e6,80e6,5);
    fprintf('原始波形已降采样 %d 倍，剩余点数: %d\n', Downsample_Factor, length(raw_sig));
end

%% === 3. 绘图 (4行 x 2列) ===
f = figure('Units', 'pixels', 'Color', 'w', 'Position', [100, 100, 1000, 900]);
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 绘图参数
pt_size = 5;        
alpha_val = 0.8;    

% 定义列标题，用于辅助说明
col_titles = {'(Upsample: 2)', '(Upsample: 50)'};

% --- 第 1 行：波形图 (Waveform) ---
for i = 1:2
    nexttile;
    if has_waveform
        plot(t_raw, raw_sig, 'LineWidth', 0.5, 'Color', 'b');
        xlim([0, T_us_max]);
        grid on;
        if i==1, ylabel('Amplitude'); end 
        title(['Waveform ' col_titles{i}]);
    else
        text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center');
    end
end

% --- 第 2 行：方位角 vs 时间 (Azimuth vs Time) ---
% Col 1
nexttile;
plot_time_scatter(data_1, 'Azimuth', pt_size, alpha_val, T_us_max);
ylabel('Azimuth (°)'); title('Azimuth vs Time');
ylim([125, 170]);

% Col 2
nexttile;
plot_time_scatter(data_2, 'Azimuth', pt_size, alpha_val, T_us_max);
title('Azimuth vs Time');
ylim([125, 170]);

% --- 第 3 行：仰角 vs 时间 (Elevation vs Time) ---
% Col 1
nexttile;
plot_time_scatter(data_1, 'Elevation', pt_size, alpha_val, T_us_max);
ylabel('Elevation (°)'); title('Elevation vs Time');
ylim([5, 45]);

% Col 2
nexttile;
plot_time_scatter(data_2, 'Elevation', pt_size, alpha_val, T_us_max);
title('Elevation vs Time');
ylim([5, 45]);

% --- 第 4 行：仰角 vs 方位角 (Elevation vs Azimuth) ---
% Col 1
nexttile;
plot_az_el(data_1, pt_size, alpha_val);
ylabel('Elevation (°)'); xlabel('Azimuth (°)'); 
title('El vs Az');
xlim([125, 170]);
ylim([5, 45]);

% Col 2
nexttile;
plot_az_el(data_2, pt_size, alpha_val);
xlabel('Azimuth (°)'); 
title('El vs Az');
xlim([125, 170]);
ylim([5, 45]);

% --- 全局设置 ---
% 统一 Colorbar
cb = colorbar;
cb.Layout.Tile = 'east';
ylabel(cb, 'Time (\mus)', 'FontSize', 11);
colormap(jet);
caxis([0, T_us_max]); 

% 统一坐标轴风格
all_axes = findall(f, 'type', 'axes');
set(all_axes, 'FontSize', 10, 'LineWidth', 1.0, 'Box', 'on', ...
    'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15]);

fprintf('绘图完成。\n点数统计:\n Upsample 2:  %d\n Upsample 50: %d\n', ...
    height(data_1), height(data_2));

%% === 辅助函数 ===
function filteredT = read_and_filter(fname, base_loc, t123, rcorr_default, T_us_max)
    if ~isfile(fname)
        warning(['文件不存在: ' fname]);
        filteredT = table();
        return;
    end
    
    opts = detectImportOptions(fname);
    opts.VariableNamingRule = 'preserve';
    T = readtable(fname, opts);
    
    % 1. 基本范围筛选
    max_samples = T_us_max * 200; 
    idx = T.Start_loc > base_loc & T.Start_loc < (base_loc + max_samples) ;
    if ismember('Azimuth', T.Properties.VariableNames) && ismember('Elevation', T.Properties.VariableNames)
        
        bad_region = (T.Azimuth > 160) & (T.Azimuth < 250) & ...
                     (T.Elevation > 0) & (T.Elevation < 50);
                 
        idx = idx & (~bad_region);
    end
    
    % 2. Rcorr 筛选
    if ismember('Rcorr', T.Properties.VariableNames)
        % 使用传入的 rcorr_default
        idx = idx & (T.Rcorr > rcorr_default);
    end

    % 3. 其他质量控制
    if ismember('t123', T.Properties.VariableNames)
        idx = idx & (abs(T.t123) < t123); 
    end
    if ismember('Elevation', T.Properties.VariableNames)
        idx = idx & (T.Elevation < 50);
    end
    if ismember('Azimuth', T.Properties.VariableNames)
        idx = idx & (T.Azimuth < 250);
    end
    
    filteredT = T(idx, :);
    filteredT.Time_us = (filteredT.Start_loc - base_loc) / 200e6 * 1e6;
end

function plot_time_scatter(data, y_field, sz, alp, t_max)
    if isempty(data), return; end
    x = data.Time_us;
    y = data.(y_field);
    scatter(x, y, sz, x, 'filled', 'MarkerFaceAlpha', alp);
    xlim([0, t_max]);
    grid on;
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);
end

function plot_az_el(data, sz, alp)
    if isempty(data), return; end
    x = data.Azimuth;
    y = data.Elevation;
    c = data.Time_us;
    scatter(x, y, sz, c, 'filled', 'MarkerFaceAlpha', alp);
    grid on;
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);
end

function signal = read_signal(signal_path, r_length, r_location)
    fid = fopen(signal_path, 'r');
    if fid == -1, error('无法打开文件'); end
    status = fseek(fid, r_location * 2, 'bof');
    if status == -1, fclose(fid); error('fseek 失败'); end
    signal = fread(fid, r_length, 'int16');
    fclose(fid);
end

function filtered_signal = filter_bp(signal,f1,f2,order)
    Fs = 200e6;
    fn = Fs/2;
    Wn = [f1 f2]/fn;
    [b,a] = butter(order,Wn); 
    filtered_signal = filtfilt(b,a,signal);
end