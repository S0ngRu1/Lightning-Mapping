clear;
clc;

% --- 配置 ---
file_path = 'E:\\LHW\\测试\\MZYP_2025-03-06_17-54-00.dat'; % 确认文件路径正确

% --- 初始化 ---
station_info = struct();
T_Year = [];
T_Month = [];
T_Day = [];
T_Hour = [];
T_Min = [];
T_Sec = [];
T_Nanosec = [];
T_PeakNum = [];

% --- 1. 打开文件 ---
fid = fopen(file_path, 'rt', 'n', 'UTF-8'); % 使用 'rt' 文本模式打开, 尝试UTF-8编码
if fid < 0
    error('无法打开文件: %s', file_path);
end

% --- 2. 读取第一行 (站点信息) ---
header_line = fgetl(fid);
if ~ischar(header_line)
    fclose(fid);
    error('文件为空或无法读取第一行。');
end

% --- 3. 解析站点信息 (保留您之前的代码) ---
tokens = regexp(header_line, 'Station Name:\s*([^\t]+)', 'tokens');
if ~isempty(tokens), station_info.station_name = strtrim(tokens{1}{1}); end
tokens = regexp(header_line, 'GPS Long:\s*([\d.]+)', 'tokens');
if ~isempty(tokens), station_info.GPS_Long = str2double(tokens{1}{1}); end
tokens = regexp(header_line, 'GPS Lat:\s*([\d.]+)', 'tokens');
if ~isempty(tokens), station_info.GPS_Lat = str2double(tokens{1}{1}); end
tokens = regexp(header_line, 'Sample Frequency:\s*(\d+)', 'tokens');
if ~isempty(tokens), station_info.Sample_Frequency = str2double(tokens{1}{1}); end
tokens = regexp(header_line, 'Point Number:\s*(\d+)', 'tokens');
if ~isempty(tokens), station_info.Point_Number = str2double(tokens{1}{1}); end
tokens = regexp(header_line, 'preTrigger Number:\s*(\d+)', 'tokens');
if ~isempty(tokens), station_info.preTrigger_Number = str2double(tokens{1}{1}); end
tokens = regexp(header_line, 'GPS Status:\s*(\w+)', 'tokens');
if ~isempty(tokens), station_info.GPS_status = strtrim(tokens{1}{1}); end

% --- 4. 逐行读取剩余内容 ---
rest_lines = {}; % 使用元胞数组存储每一行
while ~feof(fid) % 循环直到文件末尾
    line = fgetl(fid); % 读取一行
    if ischar(line) % 检查是否成功读取到一行 (fgetl在末尾返回-1)
        rest_lines{end+1} = line; % 将读取到的行添加到元胞数组
    end
end
fclose(fid); % 关闭文件

% --- 5. 合并行为一个字符串 ---
% 使用 'newline' 作为分隔符，这会自动使用适合您操作系统的换行符。
% 'fgetl' 会去掉换行符，所以我们需要加回来。
file_content_rest = strjoin(rest_lines, newline);

% --- 显示站点信息 ---
disp('提取到的站点信息:');
disp(station_info);
disp('------------------------------------');
% --- 调试信息 (可选) ---
% fprintf('读取并合并了 %d 行。\n', length(rest_lines));
% disp('合并后的内容 (前1000字符):');
% disp(file_content_rest(1:min(1000, end)));
% disp('------------------------------------');

% --- 6. 分割和处理 (与之前相同) ---
blocks = regexp(file_content_rest, '(Trigger Time:)', 'split');

fprintf('分割出了 %d 个潜在块。\n', length(blocks)-1); % 调试信息
disp('------------------------------------');


for i = 2:length(blocks)
    current_block_str = blocks{i};

    trimmed_block = strtrim(current_block_str);
    if isempty(trimmed_block)
        % fprintf('跳过块 %d (空或非Trigger Time开头)。\n', i); % 调试信息
        continue;
    end

    % fprintf('--- 正在处理块 %d ---\n', i); % 调试信息

    lines = strsplit(current_block_str, {'\n', '\r'}, 'CollapseDelimiters', false);
    lines = lines(~cellfun('isempty', lines));

    if isempty(lines)
        continue;
    end

    % --- 提取 Trigger Time 并分别解析 ---
    trigger_line = lines{1};
    time_match = regexp(trigger_line, '\s*(\d{4}-\d{2}-\d{2})H(\d{2})M(\d{2})S([\d.]+)', 'tokens');
    if isempty(time_match)
        warning('块 %d: 未找到有效的 Trigger Time: %s', i, trigger_line);
        continue;
    end

    % 提取各个部分
    date_str = time_match{1}{1};     % YYYY-MM-DD
    hour_str = time_match{1}{2};     % HH
    min_str = time_match{1}{3};      % MM
    sec_full_str = time_match{1}{4}; % SS.sssssssss

    % 解析年月日
    date_parts = sscanf(date_str, '%d-%d-%d');
    year_val = date_parts(1);
    month_val = date_parts(2);
    day_val = date_parts(3);

    % 解析时分
    hour_val = str2double(hour_str);
    min_val = str2double(min_str);

    % 解析秒和纳秒
    sec_parts = strsplit(sec_full_str, '.');
    sec_val = str2double(sec_parts{1});

    nanosec_val = 0; % 如果没有小数部分，纳秒为0
    if length(sec_parts) > 1 && ~isempty(sec_parts{2})
        frac_str = sec_parts{2};
        % 确保小数部分不超过9位 (纳秒级别)
        if length(frac_str) > 9
            frac_str = frac_str(1:9);
            warning('块 %d: 时间戳精度超过纳秒，已截断。', i);
        end
        % 将小数部分转换为数字，并乘以相应的10的幂次，使其成为纳秒
        nanosec_val = str2double(frac_str) * (10^(9 - length(frac_str)));
    end
    % --- 提取 Peak Number ---
    peak_number_line_idx = find(startsWith(lines, 'Peak Number:'), 1);
    if isempty(peak_number_line_idx)
        warning('块 %d: 未找到 Peak Number', i);
        continue;
    end
    peak_number_line = lines{peak_number_line_idx};
    peak_number_str = regexp(peak_number_line, 'Peak Number:\s*(\d+)', 'tokens');
    if isempty(peak_number_str)
        warning('块 %d: 无法解析 Peak Number: %s', i, peak_number_line);
        continue;
    end
    peak_number = str2double(peak_number_str{1}{1});

    % --- 存储到临时数组 (用于表格) ---
    T_Year(end+1) = year_val;
    T_Month(end+1) = month_val;
    T_Day(end+1) = day_val;
    T_Hour(end+1) = hour_val;
    T_Min(end+1) = min_val;
    T_Sec(end+1) = sec_val;
    T_Nanosec(end+1) = nanosec_val;
    T_PeakNum(end+1) = peak_number;
end
results_table = table(T_Year', T_Month', T_Day', T_Hour', T_Min', T_Sec', T_Nanosec', ...
            'VariableNames', {'Year', 'Month', 'Day', 'Hour', 'Min', 'Sec', 'Nanosec'});


grouping_vars = {'Year', 'Month', 'Day', 'Hour', 'Min', 'Sec'};
counts_per_second = groupcounts(results_table, grouping_vars);
counts_per_second = sortrows(counts_per_second, grouping_vars);
% 显示统计结果
disp('每秒的数据点数量统计结果:');
disp(counts_per_second);

% --- 显示或保存结果 ---
disp('------------------------------------');
disp('数据块提取完成。');