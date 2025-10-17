%% %%%tdms 文件读取 引雷点   
function signal = read_signal_chj(signal_path, r_length, r_location, skip_byte)
    % 读取交替存储的双通道信号，支持指定读取的采样点数量
    % 输入参数：
    %   signal_path  - 信号文件路径
    %   r_length     - 每个通道需要读取的采样点数量
    %   r_location   - 信号在文件中的起始偏移位置（字节）
    %   skip_byte    - 每个采样点的字节数（对于int16通常为2）
    % 输出参数：
    %   ch1_data, ch2_data - 两个通道的信号数据（长度均为r_length）
    
    fid = fopen(signal_path, 'rb');
    if fid == -1
        error('无法打开文件: %s', signal_path);
    end
    try
        % === 读取通道1数据（交替存储的奇数位） ===
        % 定位到信号起始位置
        fseek(fid, r_location, 'bof');
        % 读取r_length个int16数据，步长为2*skip_byte（跳过通道2的点）
        signal = fread(fid, r_length, 'int16', 2*skip_byte, 'b');
        fclose(fid);
    catch ME
        fclose(fid);
        error('读取信号失败: %s', ME.message);
    end
end
