clear
fid = fopen('E:\测试\特征数据\250407164545.070.bin', 'rb'); % 使用 'rb' 模式读取二进制文件
text = {}; % 初始化一个单元格数组，用于存储从文件读取的数据
% 判断文件是否成功打开
if fid > 0
% 初始化数组
all_data = [];
time_stamps = [];
% 循环读取直到文件结束
while ~feof(fid)
% 读取触发时间(8字节double，大端字节序)
trigger_time = fread(fid, 1, 'double', 'b');
if isempty(trigger_time)
break; % 文件结束
end
% 读取数据长度(4字节int32，大端字节序)
data_length = fread(fid, 1, 'int32', 'b');
% 读取数据(int16，大端字节序)
data = fread(fid, data_length, 'int16', 'b');
% 存储数据
time_stamps = [time_stamps; trigger_time];
all_data = [all_data, data];
end
% 关闭文件
fclose(fid);
% 重组数据为8通道格式（假设数据是8通道的）
if ~isempty(all_data)
num_channels = 8;
data = reshape(all_data, num_channels, []);
else
data = zeros(1, 8);
time_stamps = 0;
end
else
% 文件打开失败
data = zeros(1, 8);
time_stamps = 0;
end
% a = [61 61 61 5060 5060 5060 64
% ];
% a_1d = reshape(a, [], 1); % 将 a 转换为一维列向量
% column = 1;