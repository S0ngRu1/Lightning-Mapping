function [start_loc, azimuth, elevation, Rcorr, t123] = read_result(result_path, start_read_line, end_read_line)
    % 打开文件读取数据
    fileID = fopen(result_path, 'r');
    if fileID == -1
        error('无法打开文件！请检查文件路径');
    end

    % 读取表头行并忽略
    header = fgetl(fileID);

    % 读取文件中的所有数据
    data = textscan(fileID, '%f %f %f %f %f %f %f %f %f %f %f', 'Delimiter', '\t');
    
    % 关闭文件
    fclose(fileID);

    % 获取 Start_loc 列数据
    start_loc_full = data{1};

    % 找到大于 start_read_line 的第一个索引
    start_index = find(start_loc_full > start_read_line, 1, 'first');
    if isempty(start_index)
        error('start_read_line:未找到 Start_loc 大于 %d 的数据', start_read_line);
    end

    % 找到大于 end_read_line 的第一个索引
    end_index = find(start_loc_full > end_read_line, 1, 'first');
    if isempty(end_index)
        error('end_read_line:未找到 Start_loc 大于 %d 的数据', end_read_line);
    end

    % 从找到的位置开始提取数据
    start_loc = start_loc_full(start_index:end_index);          % Start_loc
    azimuth = data{8}(start_index:end_index);                   % Azimuth
    elevation = data{9}(start_index:end_index);                 % Elevation
    Rcorr = data{10}(start_index:end_index);                    % Rcorr
    t123 = data{11}(start_index:end_index);                     % t123
end