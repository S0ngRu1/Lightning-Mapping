clear; clc; close all;

%% === 1. 参数设置 ===
target_file = 'results\20240822165932_result_yld_window_ADAPTIVE_1e6_factor4.txt';

fs = 200e6;              % 采样率
Start_loc_Base = 3.93e8;  % 基准采样点位置 (Start Location)

% 【修改】最大时间范围改为“采样点”形式
Range_Samples = 7e6;     % 例如：1e8 个采样点 (对应 0.5秒 @ 200MHz)

% 筛选参数 
t123_th = 1;      % 闭合差阈值
rcorr_th = 0.3;   % 默认相关系数阈值

%% === 2. 数据读取与预处理 ===
% 读取结果文件 (内部逻辑已修改为使用采样点筛选)
data = read_and_filter(target_file, Start_loc_Base, t123_th, rcorr_th, Range_Samples);

%% === 3. 绘图 (仅绘制 Elevation vs Azimuth) ===
% 设置单张图的窗口大小
f = figure('Units', 'pixels', 'Position', [300, 200, 650, 550], 'Color', 'w');

pt_size = 8;        
alpha_val = 1;    

if ~isempty(data)
    x = data.Azimuth;
    y = data.Elevation;
    c = data.Time_us; % 颜色依然映射为相对时间(us)，方便阅读，但数据范围由采样点决定
    
    % 绘制散点图
    scatter(x, y, pt_size, c, 'filled', 'MarkerFaceAlpha', alpha_val);
    
    % --- 图形修饰 ---
    grid on;
    box on;
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);
    set(gca, 'FontSize', 12, 'LineWidth', 1.0);
    
    title('(d) Elevation vs Azimuth');
    ylabel('Elevation (°)'); 
    xlabel('Azimuth (°)');
    
    % 坐标轴范围设置 (根据需要取消注释)
    xlim([125, 170]); 
    % ylim([0, 90]);
    
    % --- Colorbar 设置 ---
    cb = colorbar;
    ylabel(cb, 'Time (\mus)', 'FontSize', 12);
    colormap(jet);
    
    % 计算对应的最大时间(us)用于颜色轴限制
    Max_Time_us = Range_Samples / fs * 1e6; 
    caxis([0, Max_Time_us]); 
    
else
    text(0.5, 0.5, 'No Data Found in this Range', 'HorizontalAlignment', 'center', 'FontSize', 14);
end

fprintf('绘图完成。筛选后的数据点数: %d\n', height(data));

%% === 辅助函数 ===
function filteredT = read_and_filter(fname, base_loc, t123, rcorr_default, range_samples)
    if ~isfile(fname)
        warning(['文件不存在: ' fname]);
        filteredT = table();
        return;
    end
    
    opts = detectImportOptions(fname);
    opts.VariableNamingRule = 'preserve';
    T = readtable(fname, opts);
    
    % 1. 基本范围筛选 (使用采样点 Range_Samples)
    % 筛选条件：Start_loc 在 [Base, Base + Range] 之间
    idx = T.Start_loc > base_loc & T.Start_loc < (base_loc + range_samples);
    
    % 2. 动态 Rcorr 筛选 
    if ismember('Rcorr', T.Properties.VariableNames)
        if ismember('Win_Len', T.Properties.VariableNames)
            mask_rcorr = false(height(T), 1);
            
            % 规则 1: Win_Len = 512
            idx_512 = (T.Win_Len == 512);
            mask_rcorr(idx_512) = T.Rcorr(idx_512) > 0.65;
            
            % 规则 2: Win_Len = 1024
            idx_1024 = (T.Win_Len == 1024);
            mask_rcorr(idx_1024) = T.Rcorr(idx_1024) > 0.35;
            
            % 规则 3: Win_Len = 2048
            idx_2048 = (T.Win_Len == 2048);
            mask_rcorr(idx_2048) = T.Rcorr(idx_2048) > 0.2;
            
            % 规则 4: Win_Len = 4096
            idx_4096 = (T.Win_Len == 4096);
            mask_rcorr(idx_4096) = T.Rcorr(idx_4096) > 0.2;
            
            % 其他长度
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
        idx = idx & (T.Elevation < 85);
    end
    if ismember('Azimuth', T.Properties.VariableNames)
        idx = idx & (T.Azimuth < 160);
    end
    
    filteredT = T(idx, :);
    
    % 计算相对时间 (微秒)，用于颜色映射
    filteredT.Time_us = (filteredT.Start_loc - base_loc) / 200e6 * 1e6;
end