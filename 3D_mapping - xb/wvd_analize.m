% --- 原始参数设定 ---
fs = 200e6;
signal_length = 60;
r_loction_yld = 4.698e8 + 23300 + 720;
channels = {'CH1', 'CH2', 'CH3'};
numChannels = length(channels);

% --- 循环处理每个通道 ---
for i = 1:numChannels
    % 读取并滤波信号 (与您原代码相同)
    filename = sprintf('..\\20240822165932.6610%s.dat', channels{i});
    signal = read_signal(filename, signal_length, r_loction_yld);
    signal = filter_bp(signal, 30e6, 80e6, 5);
    
    % % ----------【修改部分开始】----------
    % % 原CWT代码被替换
    % fb = cwtfilterbank(...);
    % [cfs, f] = wt(fb, signal);
    
    % 使用MATLAB内置的wvd函数计算平滑伪Wigner-Ville分布
    % [spwvd_map, f, t] = wvd(signal, fs, 'smoothedPseudo');
    % 'smoothedPseudo' 模式能更好地抑制交叉项，使时频图更清晰
    [spwvd_map, f_wvd, t_wvd] = wvd(signal, fs, 'smoothedPseudo');
    % ----------【修改部分结束】----------
    
    
    % --- 绘图部分 ---
    figure('Name', sprintf('通道 %s 分析结果', channels{i}));
    
    % 绘制时频图 (使用wvd的结果)
    subplot(2, 1, 1);
    
    % % 原CWT绘图代码
    % scalogram_data = abs(cfs);
    % imagesc((1:length(signal))/fs*1e6, f/1e6, scalogram_data);
    
    % 绘制 Smoothed Pseudo WVD 时频图
    imagesc(t_wvd * 1e6, f_wvd / 1e6, abs(spwvd_map)); % 时间单位转换为us, 频率单位转换为MHz
    axis xy;
    xlabel('时间 (\mus)');
    ylabel('频率 (MHz)');
    % 更新标题以反映新的分析方法
    title(sprintf('通道 %s 时频图 (Smoothed Pseudo WVD)', channels{i}));
    colorbar;
    
    % 绘制时域波形 (与您原代码相同)
    subplot(2, 1, 2);
    plot((1:length(signal))/fs*1e6, signal);
    xlabel('时间 (\mus)');
    ylabel('信号幅度');
    title(sprintf('通道 %s 时域波形', channels{i}));
    grid on;
end


