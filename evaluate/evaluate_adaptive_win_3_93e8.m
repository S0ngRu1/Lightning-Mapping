clear; clc; close all;

%% === 1. 参数设置 ===
% --- 文件路径 ---
file_col1 = '..\2024\2d\results\20240822165932_result_yld_3.65e8_5e8_window_512_128_阈值4倍标准差_去零飘_30_80_hann_with_error.txt';
file_col2 = '..\2024\2d\results\20240822165932_result_yld_3.65e8_5e8_window_4096_1024_阈值4倍标准差_去零飘_30_80_hann.txt';
file_col3 = '..\2024\2d\results\result_yld_3.65e8_5.6e8_window_ADAPTIVE_1e4_factor3.txt';

% 原始波形文件
raw_file_ch1 = '..\2024\20240822165932.6610CH1.dat'; 
has_waveform = isfile(raw_file_ch1); 

% --- 通用参数 ---
fs = 200e6;              % 采样率
Start_loc_Base = 3.93e8;  % 基准采样点位置
T_us_max = 30000;        % 最大显示时间 (us)

% --- 【关键】局部放大区域设置 ---
Zoom_Az_Lim = [140, 148];  % 方位角放大范围
Zoom_El_Lim = [33, 40];    % 仰角放大范围

% --- 全局视图范围 ---
Global_Az_Lim = [125, 170];
Global_El_Lim = [30, 40];

%% === 2. 数据读取与预处理 ===
% Col 1: 512窗口
data_1 = read_and_filter(file_col1, Start_loc_Base, 1, 0.6, T_us_max);
% Col 2: 4096窗口
data_2 = read_and_filter(file_col2, Start_loc_Base, 1, 0.1, T_us_max);
% Col 3: Adaptive
data_3 = read_and_filter(file_col3, Start_loc_Base, 1.0, 0.1, T_us_max);

% --- 读取原始波形 ---
if has_waveform
    read_len = ceil(T_us_max * 1e-6 * fs); 
    raw_sig = read_signal(raw_file_ch1, read_len, Start_loc_Base);
    t_full = (0:length(raw_sig)-1) / fs * 1e6; 
    
    % 降采样 (仅用于绘图，减少卡顿)
    Downsample_Factor = 50; 
    raw_sig = raw_sig(1:Downsample_Factor:end);
    t_raw = t_full(1:Downsample_Factor:end); 
    % 简单滤波以获得更清晰的波形
    raw_sig = filter_bp(raw_sig, 30e6, 80e6, 5);
else
    raw_sig = []; t_raw = [];
end

%% === 3. 绘图 (6行 x 3列 虚拟网格) ===
% 设置画布大小 (增加高度以容纳三行图)
f = figure('Units', 'pixels', 'Color', 'w', 'Position', [50, 50, 1200, 900]);

% --- 使用 6x3 网格来实现高度分配 (1:3:2) ---
t = tiledlayout(6, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% 绘图参数
pt_size_overview = 5;   % 全景图点大小     
pt_size_zoom = 10;      % 放大图点大小
alpha_val = 0.8;    

% 定义列标题
col_titles = {'(win:512)', '(win:4096)', '(Adaptive win)'};

% =========================================================================
% 第 1 行：信号波形图 (Waveform) [占用 1 行网格]
% =========================================================================
for i = 1:3
    % ax_sig = nexttile(索引)
    ax_sig = nexttile(i, [1, 1]); 
    
    if ~isempty(raw_sig)
        plot(t_raw, raw_sig, 'LineWidth', 0.5, 'Color', [0 0.4470 0.7410]); % 蓝色波形
        xlim([0, T_us_max]);
        grid on;
        set(gca, 'FontSize', 9);
        
        % 仅第一列显示 Y 轴标签
        if i == 1
            ylabel('Amp (a.u.)');
        else
            yticklabels([]);
        end
        
        % 【修改】添加横坐标刻度和标题
        xlabel('Time (\mus)', 'FontSize', 10); 
    end
    
    % 标题移动到最上方
    title(col_titles{i}, 'FontSize', 12, 'FontWeight', 'bold');
    % 添加编号 (a1, a2, a3)
    text(ax_sig, 0.02, 0.90, sprintf('(a%d) Waveform', i), 'Units', 'normalized', 'FontWeight', 'bold');
end

% =========================================================================
% 第 2 行：整体对比 (Overview) [占用 3 行网格]
% =========================================================================
% 起始索引 = 4

% --- 2.1 Col 1 Overview ---
ax1 = nexttile(4, [3, 1]); 
plot_overview(ax1, data_1, Zoom_Az_Lim, Zoom_El_Lim, Global_Az_Lim, Global_El_Lim, pt_size_overview, alpha_val);
ylabel('Elevation (°)'); 
text(ax1, 0.02, 0.95, '(b1) Overview', 'Units', 'normalized', 'FontWeight', 'bold');

% --- 2.2 Col 2 Overview ---
ax2 = nexttile(5, [3, 1]);
plot_overview(ax2, data_2, Zoom_Az_Lim, Zoom_El_Lim, Global_Az_Lim, Global_El_Lim, pt_size_overview, alpha_val);
text(ax2, 0.02, 0.95, '(b2) Overview', 'Units', 'normalized', 'FontWeight', 'bold');

% --- 2.3 Col 3 Overview ---
ax3 = nexttile(6, [3, 1]);
plot_overview(ax3, data_3, Zoom_Az_Lim, Zoom_El_Lim, Global_Az_Lim, Global_El_Lim, pt_size_overview, alpha_val);
text(ax3, 0.02, 0.95, '(b3) Overview', 'Units', 'normalized', 'FontWeight', 'bold');

% =========================================================================
% 第 3 行：局部放大对比 (Zoomed) [占用 2 行网格]
% =========================================================================
% 起始索引 = 13 (前3行共 1x3 + 3x3 = 12格占用)

% --- 3.1 Col 1 Zoom ---
ax4 = nexttile(13, [2, 1]);
plot_zoom(ax4, data_1, Zoom_Az_Lim, Zoom_El_Lim, pt_size_zoom, alpha_val);
ylabel('Elevation (°)'); xlabel('Azimuth (°)');
text(ax4, 0.02, 0.95, '(c1) Zoomed', 'Units', 'normalized', 'FontWeight', 'bold', 'Color', 'r');

% --- 3.2 Col 2 Zoom ---
ax5 = nexttile(14, [2, 1]);
plot_zoom(ax5, data_2, Zoom_Az_Lim, Zoom_El_Lim, pt_size_zoom, alpha_val);
xlabel('Azimuth (°)');
text(ax5, 0.02, 0.95, '(c2) Zoomed', 'Units', 'normalized', 'FontWeight', 'bold', 'Color', 'r');

% --- 3.3 Col 3 Zoom ---
ax6 = nexttile(15, [2, 1]);
plot_zoom(ax6, data_3, Zoom_Az_Lim, Zoom_El_Lim, pt_size_zoom, alpha_val);
xlabel('Azimuth (°)');
text(ax6, 0.02, 0.95, '(c3) Zoomed', 'Units', 'normalized', 'FontWeight', 'bold', 'Color', 'r');

% =========================================================================
% 全局设置
% =========================================================================
% 统一 Colorbar
cb = colorbar;
cb.Layout.Tile = 'east';
ylabel(cb, 'Time (\mus)', 'FontSize', 11);
colormap(jet);
caxis([0, T_us_max]); 

% 统一坐标轴风格
all_axes = findall(f, 'type', 'axes');
set(all_axes, 'FontSize', 10, 'LineWidth', 1.0, 'Box', 'on', ...
    'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15], 'TickDir', 'in');

fprintf('绘图完成。\n点数统计:\n Win512: %d\n Win4096: %d\n Adaptive: %d\n', ...
    height(data_1), height(data_2), height(data_3));

%% === 4. 辅助绘图函数 ===

% --- 绘制全景图 (带红框) ---
function plot_overview(ax, data, zoom_az, zoom_el, limit_az, limit_el, sz, alp)
    axes(ax); hold on; grid on;
    if ~isempty(data)
        scatter(data.Azimuth, data.Elevation, sz, data.Time_us, 'filled', 'MarkerFaceAlpha', alp);
    end
    % 画红框
    rectangle('Position', [zoom_az(1), zoom_el(1), diff(zoom_az), diff(zoom_el)], ...
              'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '-');
    xlim(limit_az); 
    ylim(limit_el);
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);
end

% --- 绘制局部放大图 (红色边框) ---
function plot_zoom(ax, data, zoom_az, zoom_el, sz, alp)
    axes(ax); hold on; grid on;
    if ~isempty(data)
        scatter(data.Azimuth, data.Elevation, sz, 'r', 'filled', 'MarkerFaceAlpha', alp);
    end
    xlim(zoom_az); 
    ylim(zoom_el);
    
    % 绘制红色矩形框
    rectangle('Position', [zoom_az(1), zoom_el(1), diff(zoom_az), diff(zoom_el)], ...
              'EdgeColor', 'black', 'LineWidth', 1, 'LineStyle', '-');
          
    % 设置网格线样式
    set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.5);
end

% --- 数据读取与筛选函数 (保持不变) ---
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
    
    % 2. 动态 Rcorr 筛选
    if ismember('Rcorr', T.Properties.VariableNames)
        if ismember('Win_Len', T.Properties.VariableNames)
            mask_rcorr = false(height(T), 1);
            
            idx_512 = (T.Win_Len == 512);
            mask_rcorr(idx_512) = T.Rcorr(idx_512) > 0.6;
            
            idx_1024 = (T.Win_Len == 1024);
            mask_rcorr(idx_1024) = T.Rcorr(idx_1024) > 0.4;
            
            idx_2048 = (T.Win_Len == 2048);
            mask_rcorr(idx_2048) = T.Rcorr(idx_2048) > 0.2;
            
            idx_4096 = (T.Win_Len == 4096);
            mask_rcorr(idx_4096) = T.Rcorr(idx_4096) > 0.1;
            
            other_lens = ~ismember(T.Win_Len, [512, 1024, 2048, 4096]);
            if any(other_lens)
                mask_rcorr(other_lens) = T.Rcorr(other_lens) > rcorr_default;
            end
            
            idx = idx & mask_rcorr;
        else
            idx = idx & (T.Rcorr > rcorr_default);
        end
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

% --- 文件读取辅助函数 ---
function signal = read_signal(signal_path, r_length, r_location)
    fid = fopen(signal_path, 'r');
    if fid == -1, error('无法打开文件'); end
    status = fseek(fid, r_location * 2, 'bof');
    if status == -1, fclose(fid); error('fseek 失败'); end
    signal = fread(fid, r_length, 'int16');
    fclose(fid);
end

% --- 滤波器辅助函数 ---
function filtered_signal = filter_bp(signal,f1,f2,order)
    Fs = 200e6;
    fn = Fs/2;
    Wn = [f1 f2]/fn;
    [b,a] = butter(order,Wn); 
    filtered_signal = filtfilt(b,a,signal);
end