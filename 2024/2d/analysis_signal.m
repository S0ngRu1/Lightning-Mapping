%% 正先导信号
% 1. 第一段信号（3.65 - 3.66e8 区间）
% 读取信号
signal1 = read_signal('..\\20240822165932.6610CH1.dat', 0.01e8, 3.65e8);
% 1.1 计算信号的原始均值（所有元素的均值）
mean_all1 = mean(signal1);
% 1.2 筛选出“大于原始均值”的元素
signal_above_mean1 = signal1(signal1 > mean_all1);
% 1.3 计算筛选后元素的均值（目标均值）
mean_above_mean1 = mean(signal_above_mean1);
% 1.4 计算峰值（保留原逻辑）
peak_value1 = max(signal1);

% 打印结果（可选，验证数据）
fprintf('第一段信号（3.65 - 3.66e8）：\n');
fprintf('原始均值 = %.4f\n', mean_all1);
fprintf('大于原始均值的元素的均值 = %.4f\n', mean_above_mean1);
fprintf('峰值 = %.4f\n\n', peak_value1);


% 2. 第二段信号（3.66 - 3.67e8 区间）
% 读取信号
signal2 = read_signal('..\\20240822165932.6610CH1.dat', 0.01e8, 3.66e8);
% 2.1 计算信号的原始均值（所有元素的均值）
mean_all2 = mean(signal2);
% 2.2 筛选出“大于原始均值”的元素
signal_above_mean2 = signal2(signal2 > mean_all2);
% 2.3 计算筛选后元素的均值（目标均值）
mean_above_mean2 = mean(signal_above_mean2);
% 2.4 计算峰值（保留原逻辑）
peak_value2 = max(signal2);

% 打印结果（可选，验证数据）
fprintf('第二段信号（3.66 - 3.67e8）：\n');
fprintf('原始均值 = %.4f\n', mean_all2);
fprintf('大于原始均值的元素的均值 = %.4f\n', mean_above_mean2);
fprintf('峰值 = %.4f\n', peak_value2);


% 2. 第二段信号（3.67 - 3.7e8 区间）
% 读取信号
signal2 = read_signal('..\\20240822165932.6610CH1.dat', 0.03e8, 3.67e8);
% 2.1 计算信号的原始均值（所有元素的均值）
mean_all2 = mean(signal2);
% 2.2 筛选出“大于原始均值”的元素
signal_above_mean2 = signal2(signal2 > mean_all2);
% 2.3 计算筛选后元素的均值（目标均值）
mean_above_mean2 = mean(signal_above_mean2);
% 2.4 计算峰值（保留原逻辑）
peak_value2 = max(signal2);

% 打印结果（可选，验证数据）
fprintf('第二段信号（3.67 - 3.7e8）：\n');
fprintf('原始均值 = %.4f\n', mean_all2);
fprintf('大于原始均值的元素的均值 = %.4f\n', mean_above_mean2);
fprintf('峰值 = %.4f\n', peak_value2);

% 2. 第二段信号（3.7 - 3.72e8 区间）
% 读取信号
signal2 = read_signal('..\\20240822165932.6610CH1.dat', 0.02e8, 3.7e8);
% 2.1 计算信号的原始均值（所有元素的均值）
mean_all2 = mean(signal2);
% 2.2 筛选出“大于原始均值”的元素
signal_above_mean2 = signal2(signal2 > mean_all2);
% 2.3 计算筛选后元素的均值（目标均值）
mean_above_mean2 = mean(signal_above_mean2);
% 2.4 计算峰值（保留原逻辑）
peak_value2 = max(signal2);

% 打印结果（可选，验证数据）
fprintf('第二段信号（3.7 - 3.72e8）：\n');
fprintf('原始均值 = %.4f\n', mean_all2);
fprintf('大于原始均值的元素的均值 = %.4f\n', mean_above_mean2);
fprintf('峰值 = %.4f\n', peak_value2);



%% 负先导信号
% 1. 第1段信号（3.816e8-3.9e8 区间，高度5000-6500）
% 读取信号
signal2 = read_signal('..\\20240822165932.6610CH1.dat', 0.084e8, 3.817e8);
% 2.1 计算信号的原始均值（所有元素的均值）
mean_all2 = mean(signal2);
% 2.2 筛选出“大于原始均值”的元素
signal_above_mean2 = signal2(signal2 > 280);
% 2.3 计算筛选后元素的均值（目标均值）
mean_above_mean2 = mean(signal_above_mean2);
% 2.4 计算峰值（保留原逻辑）
peak_value2 = max(signal2);

% 打印结果（可选，验证数据）
fprintf('第二段信号（3.816e8-3.9e8 区间，高度5000-6500）：\n');
fprintf('原始均值 = %.4f\n', mean_all2);
fprintf('大于原始均值的元素的均值 = %.4f\n', mean_above_mean2);
fprintf('峰值 = %.4f\n', peak_value2);



function signal = read_signal(signal_path, r_length,r_loction)
fid  = fopen(signal_path,'r');%读取数据的位置

%使用fseek函数将文件指针移动到指定位置，以便读取数据。
%这里指定移动位置为r_location，表示移动到指定位置开始读取数据。
fseek(fid,r_loction*2,'bof');
%使用fread函数从文件中读取数据，读取的数据长度为r_length，数据以int16格式读取。
%将读取到的数据分别保存到变量ch_1、ch_2和ch_3中。
signal = fread(fid,r_length,'int16');
%关闭所有文件
fclose(fid);
end

