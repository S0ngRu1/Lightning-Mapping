clear; clc; close all;

%% === 1. 参数设置 ===
% 这里默认使用您之前的 file_new (优化后的结果)
target_file = '..\2024\2d\results\20240822165932_result_yld_window_ADAPTIVE_1e6_factor4.txt';

% 原始波形文件
raw_file_ch1 = '..\2024\2024 822 85933.651462CH1.dat'; 
has_waveform = isfile(raw_file_ch1); 

fs = 200e6;              % 采样率
Start_loc_Base = 3.9e8;  % 基准采样点位置
T_us_max = 500e3;        % 最大时间范围 (微秒)
Downsample_Factor = 50;  % 降采样倍数

% 筛选参数 
% 注意：如果文件中有 Win_Len 列，rcorr_th 将作为"不在规则内的数据"的默认阈值
t123_th = 1;      % 闭合差阈值
rcorr_th = 0.3;   % 默认相关系数阈值 (用于无 Win_Len 列或未定义的 Win_Len)

%% === 2. 数据读取与预处理 ===
% 读取结果文件 (内部已包含动态筛选逻辑)
data = read_and_filter(target_file, Start_loc_Base, t123_th, rcorr_th, T_us_max);

% 读取原始波形 (如果有)
if has_waveform
    read_len = ceil(T_us_max * 1e-6 * fs); 
    raw_sig = read_signal(raw_file_ch1, read_len, Start_loc_Base);
    t_full = (0:length(raw_sig)-1) / fs * 1e6; 
    
    % 降采样
    raw_sig = raw_sig(1:Downsample_Factor:end);
    t_raw = t_full(1:Downsample_Factor:end); 
    fprintf('原始波形已降采样 %d 倍，剩余点数: %d\n', Downsample_Factor, length(raw_sig));
end

%% === 3. 绘图 (4行 x 1列) ===
f = figure('Units', 'pixels', 'Position', [150, 50, 700, 1100], 'Color', 'w');
t = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

pt_size = 4;        
alpha_val = 1;    

% --- 子图 1: 波形图 ---
nexttile; 
if has_waveform
    plot(t_raw, raw_sig, 'LineWidth', 0.6, 'Color', 'b');
    xlim([0, T_us_max]);
    title('(a) Time Waveform');
    ylabel('Amplitude');
    grid on;
else
    text(0.5, 0.5, 'Waveform Data Not Found', 'HorizontalAlignment', 'center');
end

% --- 子图 2: Azimuth vs Time ---
nexttile;
plot_time_scatter(data, 'Azimuth', 8, alpha_val, T_us_max);
title('(b) Azimuth vs Time');
ylabel('Azimuth (°)');

% --- 子图 3: Elevation vs Time ---
nexttile;
plot_time_scatter(data, 'Elevation', 8, alpha_val, T_us_max);
title('(c) Elevation vs Time');
ylabel('Elevation (°)');

% --- 子图 4: Elevation vs Azimuth ---
nexttile;
plot_az_el(data, pt_size, alpha_val);
title('(d) Elevation vs Azimuth');
ylabel('Elevation (°)'); xlabel('Azimuth (°)');
xlim([125, 250]); 

% --- 全局设置 ---
cb = colorbar;
cb.Layout.Tile = 'east';
ylabel(cb, 'Time (\mus)', 'FontSize', 11);
colormap(jet);
caxis([0, T_us_max]); 

all_axes = findall(f, 'type', 'axes');
set(all_axes, 'FontSize', 11, 'LineWidth', 1.0, 'Box', 'on', ...
    'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15]);

fprintf('绘图完成。筛选后的数据点数: %d\n', height(data));

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
    idx = T.Start_loc > base_loc & T.Start_loc < (base_loc + max_samples);
    
    % 2. 动态 Rcorr 筛选 
    if ismember('Rcorr', T.Properties.VariableNames)
        if ismember('Win_Len', T.Properties.VariableNames)
            % 创建一个全 false 的掩码
            mask_rcorr = false(height(T), 1);
            
            % 规则 1: Win_Len = 512, Rcorr > 0.7
            idx_512 = (T.Win_Len == 512);
            mask_rcorr(idx_512) = T.Rcorr(idx_512) > 0.7;
            
            % 规则 2: Win_Len = 1024, Rcorr > 0.6
            idx_1024 = (T.Win_Len == 1024);
            mask_rcorr(idx_1024) = T.Rcorr(idx_1024) > 0.4;
            
            % 规则 3: Win_Len = 2048, Rcorr > 0.5
            idx_2048 = (T.Win_Len == 2048);
            mask_rcorr(idx_2048) = T.Rcorr(idx_2048) > 0.2;
            
            % 规则 4: Win_Len = 4096, Rcorr > 0.3
            idx_4096 = (T.Win_Len == 4096);
            mask_rcorr(idx_4096) = T.Rcorr(idx_4096) > 0.2;
            
            % 处理不在上述规则中的 Win_Len (如果有其他长度，应用默认 rcorr)
            other_lens = ~ismember(T.Win_Len, [512, 1024, 2048, 4096]);
            if any(other_lens)
                mask_rcorr(other_lens) = T.Rcorr(other_lens) > rcorr_default;
                fprintf('注意: 发现 %d 个点的 Win_Len 不在规则(512/1024/2048/4096)中，已使用默认阈值 %.2f\n', ...
                    sum(other_lens), rcorr_default);
            end
            
            % 应用 Rcorr 筛选
            idx = idx & mask_rcorr;
        else
            % 如果没有 Win_Len 列，使用固定阈值
            idx = idx & (T.Rcorr > rcorr_default);
        end
    end

    % 3. 其他质量控制
    if ismember('t123', T.Properties.VariableNames)
        idx = idx & (abs(T.t123) < t123); 
    end
    if ismember('Elevation', T.Properties.VariableNames)
        idx = idx & (T.Elevation < 85);
    end
    if ismember('Azimuth', T.Properties.VariableNames)
        idx = idx & (T.Azimuth < 250);
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
    c = data.Time_us;
    scatter(x, y, sz, c, 'filled', 'MarkerFaceAlpha', alp);
    grid on;
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);
    
    if ~isempty(x)
        % xlim([min(x)-5, max(x)+5]);
        % ylim([min(y)-5, max(y)+5]);
    end
end

function signal = read_signal(signal_path, r_length, r_location)
    fid = fopen(signal_path, 'r');
    if fid == -1, error('无法打开文件'); end
    status = fseek(fid, r_location * 2, 'bof');
    if status == -1, fclose(fid); error('fseek 失败'); end
    signal = fread(fid, r_length, 'int16');
    fclose(fid);
end