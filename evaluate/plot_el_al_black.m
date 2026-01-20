clear; clc; close all;

%% === 1. 参数设置 ===
% --- 文件路径 ---
% 1. 传统方法 (512点窗口)
file_col1 = '..\2024\2d\results\20240822165932_result_yld_3.65e8_5e8_window_512_128_阈值4倍标准差_去零飘_30_80_hann_with_error.txt';
% 2. 新方法 (4096点窗口)
file_col2 = '..\2024\2d\results\20240822165932_result_yld_3.65e8_5e8_window_4096_1024_阈值4倍标准差_去零飘_30_80_hann.txt';
% 3. 自适应方法 (Adaptive / Loop)
file_col3 = '..\2024\2d\results\result_yld_3.65e8_5.6e8_window_ADAPTIVE_1e4_factor3.txt';

% --- 通用参数 ---
fs = 200e6;              % 采样率
Start_loc_Base = 3.93e8;  % 基准采样点位置
T_us_max = 35000;        % 最大显示时间 (us)

%% === 2. 数据读取与预处理 ===
% read_and_filter(文件名, 基准点, 闭合差阈值, Rcorr阈值, 最大时间)
% Col 1: 512窗口
data_1 = read_and_filter(file_col1, Start_loc_Base, 1, 0.6, T_us_max);
% Col 2: 4096窗口
data_2 = read_and_filter(file_col2, Start_loc_Base, 1, 0.1, T_us_max);
% Col 3: Adaptive
data_3 = read_and_filter(file_col3, Start_loc_Base, 1.0, 0.1, T_us_max);

%% === 3. 绘图 (1行 x 3列) ===
% 调整窗口形状为长条形
f = figure('Units', 'pixels', 'Position', [100, 200, 1200, 400], 'Color', 'w');
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% 绘图参数
pt_size = 5;        
alpha_val = 0.8;    % 透明度 (0-1)，0.8为深黑，越低越淡
col_titles = {'(win:512)', '(win:4096)', '(Adaptive win)'};

% --- Col 1: 传统方法 ---
nexttile;
plot_az_el_black(data_1, pt_size, alpha_val);
ylabel('Elevation (°)'); xlabel('Azimuth (°)'); 
title(['El vs Az ' col_titles{1}]);
xlim([140, 148]);
ylim([30, 42]);

% --- Col 2: 新方法 ---
nexttile;
plot_az_el_black(data_2, pt_size, alpha_val);
xlabel('Azimuth (°)'); 
title(['El vs Az ' col_titles{2}]);
xlim([140, 148]);
ylim([30, 42]);

% --- Col 3: 自适应方法 ---
nexttile;
plot_az_el_black(data_3, pt_size, alpha_val);
xlabel('Azimuth (°)'); 
title(['El vs Az ' col_titles{3}]);
xlim([140, 148]);
ylim([30, 42]);

% --- 全局设置 ---
% 统一坐标轴风格
all_axes = findall(f, 'type', 'axes');
set(all_axes, 'FontSize', 12, 'LineWidth', 1.0, 'Box', 'on', ...
    'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15]);

fprintf('绘图完成。\n点数统计: Std=%d, New=%d, Adapt=%d\n', ...
    height(data_1), height(data_2), height(data_3));

%% === 辅助函数 ===

% 【修改】绘制纯黑色散点图的函数
function plot_az_el_black(data, sz, alp)
    if isempty(data)
        text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center');
        return; 
    end
    x = data.Azimuth;
    y = data.Elevation;
    
    % 使用 'k' 代表黑色，移除颜色映射
    scatter(x, y, sz, 'red', 'filled', 'MarkerFaceAlpha', alp);
    
    grid on;
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);
end

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
            % === 动态逻辑 ===
            mask_rcorr = false(height(T), 1);
            
            % Win=512 -> >0.6
            idx_512 = (T.Win_Len == 512);
            mask_rcorr(idx_512) = T.Rcorr(idx_512) > 0.6;
            
            % Win=1024 -> >0.4
            idx_1024 = (T.Win_Len == 1024);
            mask_rcorr(idx_1024) = T.Rcorr(idx_1024) > 0.45;
            
            % Win=2048 -> >0.2
            idx_2048 = (T.Win_Len == 2048);
            mask_rcorr(idx_2048) = T.Rcorr(idx_2048) > 0.2;
            
            % Win=4096 -> >0.1
            idx_4096 = (T.Win_Len == 4096);
            mask_rcorr(idx_4096) = T.Rcorr(idx_4096) > 0.1;
            
            % 其他 -> 默认
            other_lens = ~ismember(T.Win_Len, [512, 1024, 2048, 4096]);
            if any(other_lens)
                mask_rcorr(other_lens) = T.Rcorr(other_lens) > rcorr_default;
            end
            idx = idx & mask_rcorr;
        else
            % === 固定逻辑 ===
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