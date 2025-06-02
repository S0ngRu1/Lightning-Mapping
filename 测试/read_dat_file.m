function [timestamps, frameData] = readDataFile(filename)
    % 打开文件
    fid = fopen(filename, 'r');
    
    % 检查文件是否成功打开
    if fid == -1
        error('无法打开文件');
    end
    
    % 初始化数组以存储数据
    timestamps = [];
    frameData = [];
    
    try
        % 读取文件直到结束
        while ~feof(fid)
            % 读取时间戳 (6个U32)
            timestamp = fread(fid, 6, 'uint32');
            
            % 如果读取的时间戳不完整，说明到达文件末尾
            if length(timestamp) < 6
                break;
            end
            
            % 读取frame数据 (204800个I16)
            frame = fread(fid, 10240, 'int16');
            
            % 如果读取的frame数据不完整，说明到达文件末尾
            if length(frame) < 10240
                break;
            end
            
            % 将数据添加到结果数组中
            timestamps = [timestamps; timestamp'];
            frameData = [frameData; frame'];
        end
        
    catch ME
        % 关闭文件
        fclose(fid);
        rethrow(ME);
    end
    
    % 关闭文件
    fclose(fid);
    
    % 将时间戳转换为更易读的格式
    time_struct = struct('ns', timestamps(:,1), ...
                        'us', timestamps(:,2), ...
                        'ms', timestamps(:,3), ...
                        's', timestamps(:,4), ...
                        'min', timestamps(:,5), ...
                        'hour', timestamps(:,6));
end