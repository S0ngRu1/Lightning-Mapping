% 读取数据文件
[timestamps, frameData] = read_dat_file('two_hour.dat');

timeInNanoseconds = timestamps(:,1) + ...          % 纳秒
                    timestamps(:,2) * 1e3 + ...     % 微秒
                    timestamps(:,3) * 1e6 + ...     % 毫秒
                    timestamps(:,4) * 1e9 + ...     % 秒
                    timestamps(:,5) * 60 * 1e9 + ...% 分
                    timestamps(:,6) * 3600 * 1e9;   % 时

% 计算相邻时间的差值
timeDifferences = diff(timeInNanoseconds);
timeDifferences = timeDifferences-1e9;

% 初始化结果数组
resultValues = zeros(size(timestamps(:,1)));
resultValues(1) = timestamps(1,1);
for i = 1:length(timeDifferences)
    resultValues(i+1) = resultValues(i) + timeDifferences(i);
end
plot(resultValues)

% % 绘制时间差值图
% figure;
% bar(abs(timeDifferences), 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'k');
% xlabel('样本索引');
% ylabel('时间差 (纳秒)');
% title('相邻时间差值图');
% grid on;

% timestamps 是一个 N×6 的矩阵，每行包含一个时间戳的6个分量 [ns, us, ms, s, min, hour]
% frameData 是一个 N×204800 的矩阵，每行包含一个frame的数据
% plot(frameData(1,:));
% title('First Frame Waveform');
% xlabel('Sample');
% ylabel('Amplitude');


% % 将所有frame数据连接成一个长向量
% allData = reshape(frameData', [], 1);
% % 绘制完整波形
% plot(allData);
% title('Complete Waveform');
% xlabel('Sample');
% ylabel('Amplitude');
% 
% % 添加网格便于观察
% grid on;