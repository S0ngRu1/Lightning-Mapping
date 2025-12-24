clear; clc; close all;

%% === 1. 参数设置 ===
% 文件路径
file_std = '..\2023\results\standard_2023_result.txt'; % 传统方法 (左列)
file_new = '..\2023\results\20230718175104_result_yld_3e8_6e8_window_512_128_阈值4倍标准差_去零飘_20_80_hann.txt'; % 本文方法 (右列)

% 原始波形文件 (用于绘制 a 和 e 子图)
% 如果您没有这些文件，可以将 has_waveform 置为 false
raw_file_ch1 = '..\2023\20230718175104.9180CH1.dat'; 
has_waveform = isfile(raw_file_ch1); 

fs = 200e6;              % 采样率
Start_loc_Base = 5e8;    % 基准采样点位置 (用于时间归零)
T_us_max = 400000;          % 绘图显示的最大时间范围 (微秒), 可按需调整

%% === 2. 数据读取与预处理 ===
% 读取两个结果文件
data_std = read_and_filter(file_std, Start_loc_Base, 0.001, 0.3);
data_new = read_and_filter(file_new, Start_loc_Base, 0.5, 0.6);

% 读取原始波形 (如果有)
if has_waveform
    % 读取一段数据用于展示，长度根据 T_us_max 估算
    read_len = ceil(T_us_max * 1e-6 * fs * 1.5); 
    raw_sig = read_signal(raw_file_ch1, read_len, Start_loc_Base);
    t_raw = (0:length(raw_sig)-1) / fs * 1e6; % 微秒
end

%% === 3. 绘图 (4行 x 2列) ===
f = figure('Units', 'pixels', 'Position', [100, 100, 1000, 1000], 'Color', 'w');
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 定义公共绘图参数
pt_size = 4;        % 散点大小
alpha_val = 0.6;    % 透明度

% --- 第一行：波形图 (a) & (e) ---
% (a) Raw Waveform 
nexttile; 
if has_waveform
    plot(t_raw, raw_sig, 'k', 'LineWidth', 0.5, 'Color',"blue");
    xlim([0, T_us_max]);
    title('(a) Time Waveform ');
    ylabel('Amplitude');
    grid on;
else
    text(0.5, 0.5, 'Waveform Data Not Found', 'HorizontalAlignment', 'center');
end

% (e) Waveform (Beam Steering - Ideally aligned, here showing raw for reference)
nexttile;
if has_waveform
    plot(t_raw, raw_sig, 'k', 'LineWidth', 0.5 , 'Color',"blue"); % 这里暂用原始波形示意
    xlim([0, T_us_max]);
    title('(e) Time Waveform ');
    grid on;
else
    text(0.5, 0.5, 'Waveform Data Not Found', 'HorizontalAlignment', 'center');
end

% --- 第二行：Azimuth vs Time (b) & (f) ---
% (b) Std
nexttile;
plot_time_scatter(data_std, 'Azimuth', pt_size, alpha_val, T_us_max);
title('(b) Azimuth vs Time ');
ylabel('Azimuth (°)');

% (f) New
nexttile;
plot_time_scatter(data_new, 'Azimuth', pt_size, alpha_val, T_us_max);
title('(f) Azimuth vs Time ');

% --- 第三行：Elevation vs Time (c) & (g) ---
% (c) Std
nexttile;
plot_time_scatter(data_std, 'Elevation', pt_size, alpha_val, T_us_max);
title('(c) Elevation vs Time ');
ylabel('Elevation (°)');

% (g) New
nexttile;
plot_time_scatter(data_new, 'Elevation', pt_size, alpha_val, T_us_max);
title('(g) Elevation vs Time ');

% --- 第四行：Azimuth vs Elevation (d) & (h) ---
% (d) Std - AZ-EL Plot
nexttile;
plot_az_el(data_std, pt_size, alpha_val);
title('(d) El vs Az ');
ylabel('Elevation (°)'); xlabel('Azimuth (°)');

% (h) New - AZ-EL Plot
nexttile;
plot_az_el(data_new, pt_size, alpha_val);
title('(h) El vs Az ');
xlabel('Azimuth (°)');

% --- 全局设置 ---
% 统一添加 Colorbar (放在最右侧)
cb = colorbar;
cb.Layout.Tile = 'east';
ylabel(cb, 'Time (\mus)', 'FontSize', 11);
colormap(jet);
caxis([0, T_us_max]); % 统一颜色轴范围

% 统一设置所有坐标轴风格
all_axes = findall(f, 'type', 'axes');
set(all_axes, ...
    'FontSize', 10, ...
    'LineWidth', 1.0, ...
    'Box', 'on', ...
    'XColor', [0.2 0.2 0.2], ...
    'YColor', [0.2 0.2 0.2]);


%% === 辅助函数 ===

function filteredT = read_and_filter(fname, base_loc, t123,rcorr)
    if ~isfile(fname)
        warning(['文件不存在: ' fname]);
        filteredT = table();
        return;
    end
    
    opts = detectImportOptions(fname);
    opts.VariableNamingRule = 'preserve'; % 防止变量名被修改
    T = readtable(fname, opts);
    
    % 简单的列名修正 (以防不同文件列名略有差异)
    if ismember('t12', T.Properties.VariableNames) && ~ismember('t123', T.Properties.VariableNames)
         % 假设闭合差在某些版本叫 t12... 视具体情况调整，这里假设 t123 存在
    end
    
    % 筛选逻辑
    % 1. 基本范围筛选
    idx = T.Start_loc > base_loc & T.Start_loc < (base_loc + 1e8);
    
    % 2. 质量控制 (根据您的代码)
    if ismember('t123', T.Properties.VariableNames)
        idx = idx & (abs(T.t123) < t123); % 放宽一点闭合差以便能看到数据，原代码是0.001 (1ns)
    end
    if ismember('Rcorr', T.Properties.VariableNames)
        idx = idx & (T.Rcorr > rcorr); % 原代码 0.3
    end
    if ismember('Elevation', T.Properties.VariableNames)
        idx = idx & (T.Elevation > 25 & T.Elevation < 85);
    end
    if ismember('Azimuth', T.Properties.VariableNames)
        idx = idx & (T.Azimuth > 100 & T.Azimuth < 325);
    end
    
    filteredT = T(idx, :);
    
    % 计算相对时间 (微秒)
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
    c = data.Time_us; % 颜色代表时间
    
    scatter(x, y, sz, c, 'filled', 'MarkerFaceAlpha', alp);
    grid on;
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);
    
    % 自动调整范围，稍微留点边距
    if ~isempty(x)
        margin = 5;
        xlim([min(x)-margin, max(x)+margin]);
        ylim([min(y)-margin, max(y)+margin]);
    end
end

function signal = read_signal(signal_path, r_length, r_location)
    fid = fopen(signal_path, 'r');
    if fid == -1
        error('无法打开文件');
    end
    % 移动指针到指定位置 (每个采样点2字节 int16)
    status = fseek(fid, r_location * 2, 'bof');
    if status == -1
        fclose(fid);
        error('fseek 失败');
    end
    signal = fread(fid, r_length, 'int16');
    fclose(fid);
end