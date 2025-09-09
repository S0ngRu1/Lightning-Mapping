% 设定文件路径
dataFile = 'E:\\2023data\\SignalNoise\\pulsematchingsignalnoise\\signal\\4-6e8_R0.7_t1.mat';  % 数据文件的路径
signalFile = 'E:\\2023data\\20230718175104.9180CH1.dat';  % 信号文件的路径
outputFolder = 'E:\\2023data\\datasets\\signal';  % 输出文件夹的路径

% 读取数据文件（.mat文件）
mat_data = load(dataFile);  % 加载 .mat 文件
start = mat_data.filteredTable1;
% var_names = fieldnames(start);
% first_column = start(:, 1); 
startPositions = start(:, 1);  
signalStartPositions = table2array(startPositions);

% 假设 signalStartPositions 是一个包含信号起始位置的向量
% 假设 readSignal 是一个读取信号的函数，接受起始位置和长度作为参数
% 并且返回读取的信号数据
segmentLength = 1024; % 每段信号的长度


% 确保输出文件夹存在
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% 遍历每个信号起始位置
for i = 1:length(signalStartPositions)
    startPos = signalStartPositions(i); % 获取当前的信号起始位置
    
    % 读取信号段
    signalSegment = read_signal(signalFile, segmentLength,startPos);
    
    % 生成输出文件的名称
    outputFileName = fullfile(outputFolder, sprintf('signal_segment_%d.mat', i));
    
    % 将读取到的信号数据保存到 .mat 文件中
    save(outputFileName, 'signalSegment');
end

disp('信号段保存完成。');

