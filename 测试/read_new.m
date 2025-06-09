clear;
clc;

% --- 配置 ---
file_path = 'E:\测试\2025-5-30\QKY_2025-05-30_18-03-00.bin'; % 文件路径
read_limit = 1000; % 设置读取的数据包数量上限，用于测试
num_channels = 8; % 假设数据通道数，用于重塑数据

% --- 主程序 ---
% 判断文件是否存在
if ~exist(file_path, 'file')
    error('文件不存在: %s', file_path);
end

% 使用 'rb' 模式打开二进制文件进行读取 (r=read, b=binary)
fid = fopen(file_path, 'rb');

% 判断文件是否成功打开
if fid < 0
    error('无法打开文件: %s', file_path);
else
    fprintf('文件已成功打开: %s\n', file_path);
end

% --- 初始化 ---
% 使用结构体数组来存储每个数据包的完整信息，更清晰
all_packets = struct('station_name', {}, 'longitude', {}, 'latitude', {}, ...
                     'sample_rate', {}, 'num_points', {}, 'pre_trigger', {}, ...
                     'gps_status', {}, 'data', {});
count = 0; % 数据包计数器

% --- 循环读取文件 ---
while ~feof(fid)
    % 检查是否达到读取上限
    if count >= read_limit
        fprintf('已达到读取上限 %d 个数据包。\n', read_limit);
        break;
    end
    
    % --- 读取数据包头部 ---
    
    % 1. 读取站名长度 (4字节, 无符号32位, 大端)
    name_len = fread(fid, 1, 'uint32', 'b');
    if isempty(name_len)
        disp('已到达文件末尾。');
        break; % 如果读不到长度，说明文件结束了
    end
    
    % 2. 读取站名 (char*1，长度为 name_len)
    % 'char*1=>char' 将ASCII码转换为字符，最后的 ' 转置为行向量
    station_name = fread(fid, name_len, 'char*1=>char')';
    
    % 3. 读取经度 (8字节, double, 大端)
    longitude = fread(fid, 1, 'double', 'b');
    
    % 4. 读取纬度 (8字节, double, 大端)
    latitude = fread(fid, 1, 'double', 'b');
    
    % 5. 读取采样率 (4字节, uint32, 大端)
    sample_rate = fread(fid, 1, 'uint32', 'b');
    
    % 6. 读取采样点数 (4字节, uint32, 大端)
    num_points = fread(fid, 1, 'uint32', 'b');
    
    % 7. 读取预触发点数 (4字节, uint32, 大端)
    pre_trigger = fread(fid, 1, 'uint32', 'b');
    
    % 8. 读取GPS状态 (1字节, uint8)
    gps_status = fread(fid, 1, 'uint8', 'b');
    
    % --- 读取数据内容 ---
    
    % 9. 读取数据 (int16, 大端)
    % 总数据量 = 点数 * 通道数
    total_values_to_read = num_points * num_channels;
    data_raw = fread(fid, total_values_to_read, 'int16', 'b');
    
    % 检查是否成功读取了完整的数据包
    if numel(data_raw) < total_values_to_read
        fprintf('警告: 文件末尾数据不完整，已读取 %d 个值，期望 %d 个。\n', numel(data_raw), total_values_to_read);
        break;
    end
    
    % 将一维数据重塑为多通道格式 (通道数 x 点数)
    data_reshaped = reshape(data_raw, num_channels, num_points);
    
    % --- 存储当前数据包 ---
    current_packet.station_name = station_name;
    current_packet.longitude = longitude;
    current_packet.latitude = latitude;
    current_packet.sample_rate = sample_rate;
    current_packet.num_points = num_points;
    current_packet.pre_trigger = pre_trigger;
    current_packet.gps_status = gps_status;
    current_packet.data = data_reshaped;
    
    % 将当前数据包追加到结构体数组中
    all_packets = [all_packets; current_packet];
    
    count = count + 1;
end

% --- 清理工作 ---
fclose(fid);
fprintf('文件已关闭。总共读取了 %d 个完整的数据包。\n', count);

% --- 使用示例 ---
if ~isempty(all_packets)
    % 访问第一个数据包的信息
    first_packet = all_packets(1);
    fprintf('\n--- 第一个数据包示例 ---\n');
    fprintf('站名: %s\n', first_packet.station_name);
    fprintf('经度: %f\n', first_packet.longitude);
    fprintf('纬度: %f\n', first_packet.latitude);
    fprintf('采样率: %d Hz\n', first_packet.sample_rate);
    fprintf('采样点数: %d\n', first_packet.num_points);
    fprintf('GPS状态 (1=锁定): %d\n', first_packet.gps_status);
    
    % 获取第一个数据包的第1个通道的数据
    channel_1_data = first_packet.data(1, :);
    
    % 绘制第一个数据包的第1个通道波形
    figure;
    plot(channel_1_data);
    title(sprintf('数据包 1, 通道 1 波形 (站名: %s)', first_packet.station_name));
    xlabel('采样点');
    ylabel('幅值');
    grid on;
    
    % 您也可以轻松提取所有数据包的某个字段
    all_longitudes = [all_packets.longitude];
else
    disp('未能从文件中读取任何完整的数据包。');
end