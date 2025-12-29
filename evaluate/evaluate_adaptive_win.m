clear; clc; close all;

%% === 1. 参数设置 ===
% --- 文件路径 ---
% 1. 传统方法 (512点窗口)
file_col1 = '..\2024\2d\results\20240822165932_result_yld_3.65e8_5e8_window_512_128_阈值4倍标准差_去零飘_30_80_hann_with_error.txt';
% 2. 新方法 (4096点窗口)
file_col2 = '..\2024\2d\results\20240822165932_result_yld_3.65e8_5e8_window_4096_1024_阈值4倍标准差_去零飘_30_80_hann.txt';
% 3. 自适应方法 (Adaptive / Loop)
file_col3 = '..\2024\2d\results\20240822165932_result_yld_window_ADAPTIVE_1e6_factor4.txt';

% 原始波形文件
raw_file_ch1 = '..\2024\2024 822 85933.651462CH1.dat'; 
has_waveform = isfile(raw_file_ch1); 

% --- 通用参数 ---
fs = 200e6;              % 采样率
Start_loc_Base = 3.9e8+76000000;  % 基准采样点位置
T_us_max = 30e3;        % 最大显示时间 (us)
Downsample_Factor = 50;  % 波形降采样倍数

%% === 2. 数据读取与预处理 ===
% 分别读取三列数据，应用各自的筛选阈值
% read_and_filter(文件名, 基准点, 闭合差阈值, Rcorr阈值, 最大时间)

% Col 1: 512窗口 (t123=0.5, Rcorr=0.7)
data_1 = read_and_filter(file_col1, Start_loc_Base, 0.5, 0.7, T_us_max);

% Col 2: 4096窗口 (t123=1.0, Rcorr=0.2)
data_2 = read_and_filter(file_col2, Start_loc_Base, 1.0, 0.3, T_us_max);

% Col 3: Adaptive (t123=1.0, Rcorr=0.3默认) 
% *注：函数内部会检测 Win_Len 列，若存在则使用动态阈值，Rcorr=0.3 仅作为备用
data_3 = read_and_filter(file_col3, Start_loc_Base, 1.0, 0.1, T_us_max);

% --- 读取原始波形 ---
if has_waveform
    read_len = ceil(T_us_max * 1e-6 * fs); 
    raw_sig = read_signal(raw_file_ch1, read_len, Start_loc_Base);
    t_full = (0:length(raw_sig)-1) / fs * 1e6; 
    
    % 降采样
    raw_sig = raw_sig(1:Downsample_Factor:end);
    t_raw = t_full(1:Downsample_Factor:end); 
    fprintf('原始波形已降采样 %d 倍，剩余点数: %d\n', Downsample_Factor, length(raw_sig));
end

%% === 3. 绘图 (4行 x 3列) ===
f = figure('Units', 'pixels', 'Position', [50, 50, 1400, 1000], 'Color', 'w');
t = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% 绘图参数
pt_size = 2;        
alpha_val = 0.8;    

% 定义列标题，用于辅助说明
col_titles = {'(win:512)', '(win:4096)', '(Adaptive win)'};

% --- 第 1 行：波形图 (Waveform) ---
% 为了布局对齐，我们在三列都绘制波形 (或者您可以只画中间)
for i = 1:3
    nexttile;
    if has_waveform
        plot(t_raw, raw_sig, 'LineWidth', 0.5, 'Color', 'b');
        xlim([0, T_us_max]);
        grid on;
        if i==1, ylabel('Amplitude'); end % 只在第一列显示Y轴标签
        title(['Waveform ' col_titles{i}]);
    else
        text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center');
    end
end

% --- 第 2 行：方位角 vs 时间 (Azimuth vs Time) ---
% Col 1
nexttile;
plot_time_scatter(data_1, 'Azimuth', pt_size, alpha_val, T_us_max);
ylabel('Azimuth (°)'); title('Azimuth vs Time ');

% Col 2
nexttile;
plot_time_scatter(data_2, 'Azimuth', pt_size, alpha_val, T_us_max);
title('Azimuth vs Time ');

% Col 3
nexttile;
plot_time_scatter(data_3, 'Azimuth', pt_size, alpha_val, T_us_max);
title('Azimuth vs Time ');

% --- 第 3 行：仰角 vs 时间 (Elevation vs Time) ---
% Col 1
nexttile;
plot_time_scatter(data_1, 'Elevation', pt_size, alpha_val, T_us_max);
ylabel('Elevation (°)'); title('Elevation vs Time ');
ylim([0, 82]);

% Col 2
nexttile;
plot_time_scatter(data_2, 'Elevation', pt_size, alpha_val, T_us_max);
title('Elevation vs Time');
ylim([0, 82]);

% Col 3
nexttile;
plot_time_scatter(data_3, 'Elevation', pt_size, alpha_val, T_us_max);
title('Elevation vs Time');
ylim([0, 82]);

% --- 第 4 行：仰角 vs 方位角 (Elevation vs Azimuth) ---
% Col 1
nexttile;
plot_az_el(data_1, pt_size, alpha_val);
ylabel('Elevation (°)'); xlabel('Azimuth (°)'); 
title('El vs Az');
xlim([125, 250]);
ylim([0, 82]);
% Col 2
nexttile;
plot_az_el(data_2, pt_size, alpha_val);
xlabel('Azimuth (°)'); 
title('El vs Az ');
xlim([125, 250]);
ylim([0, 82]);
% Col 3
nexttile;
plot_az_el(data_3, pt_size, alpha_val);
xlabel('Azimuth (°)'); 
title('El vs Az ');
xlim([125, 250]);
ylim([0, 82]);
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

fprintf('绘图完成。\n点数统计: Std=%d, New=%d, Adapt=%d\n', ...
    height(data_1), height(data_2), height(data_3));

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
    % T_us_max 是微秒, fs=200MHz -> 200 samples/us
    % 稍微放宽一点上限以防边界数据丢失
    max_samples = T_us_max * 200 * 1.05; 
    idx = T.Start_loc > base_loc & T.Start_loc < (base_loc + max_samples) ;
    if ismember('Azimuth', T.Properties.VariableNames) && ismember('Elevation', T.Properties.VariableNames)
        
        bad_region = (T.Azimuth > 160) & (T.Azimuth < 250) & ...
                     (T.Elevation > 0) & (T.Elevation < 20);
                 
        % 逻辑：保留原来的 idx，并且 (AND) 还要满足“不在坏区域内” (~bad_region)
        idx = idx & (~bad_region);
    end
    % 2. 动态 Rcorr 筛选 (兼容普通文件和Adaptive文件)
    if ismember('Rcorr', T.Properties.VariableNames)
        if ismember('Win_Len', T.Properties.VariableNames)
            % === 动态逻辑 (针对 Adaptive 文件) ===
            mask_rcorr = false(height(T), 1);
            
            % Win=512 -> >0.7
            idx_512 = (T.Win_Len == 512);
            mask_rcorr(idx_512) = T.Rcorr(idx_512) > 0.55;
            
            % Win=1024 -> >0.4 (根据您之前的设定)
            idx_1024 = (T.Win_Len == 1024);
            mask_rcorr(idx_1024) = T.Rcorr(idx_1024) > 0.3;
            
            % Win=2048 -> >0.2
            idx_2048 = (T.Win_Len == 2048);
            mask_rcorr(idx_2048) = T.Rcorr(idx_2048) > 0.2;
            
            % Win=4096 -> >0.2
            idx_4096 = (T.Win_Len == 4096);
            mask_rcorr(idx_4096) = T.Rcorr(idx_4096) > 0.2;
            
            % 其他 -> 默认
            other_lens = ~ismember(T.Win_Len, [512, 1024, 2048, 4096]);
            if any(other_lens)
                mask_rcorr(other_lens) = T.Rcorr(other_lens) > rcorr_default;
            end
            
            idx = idx & mask_rcorr;
        else
            % === 固定逻辑 (针对 Std/New 文件) ===
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