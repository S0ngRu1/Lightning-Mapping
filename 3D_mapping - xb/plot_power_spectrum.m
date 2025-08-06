function plot_power_spectrum(signal, varargin)
% PLOT_POWER_SPECTRUM 绘制信号的功率谱
%   plot_power_spectrum(signal) 绘制信号的功率谱，使用默认参数
%   plot_power_spectrum(signal, Fs) 指定采样频率(Hz)
%   plot_power_spectrum(signal, Fs, title_str) 同时指定标题
%
%   输入参数：
%       signal      - 输入信号（向量）
%       Fs          - 采样频率(Hz)，默认200e6
%       title_str   - 图表标题，默认'信号功率谱'
%

    % 设置默认参数
    Fs = 200e6;             % 默认采样频率200MHz
    title_str = '信号功率谱'; % 默认标题
    
    % 解析可变参数
    if nargin >= 2
        Fs = varargin{1};
    end
    if nargin >= 3
        title_str = varargin{2};
    end
    
    % 验证输入信号
    if ~isvector(signal)
        error('输入必须是一维信号向量');
    end
    N = length(signal);
    
    % 计算功率谱（使用pwelch方法，适合长信号）
    % 自动选择合适的窗函数和分段长度
    [Pxx, f] = pwelch(signal, [], [], [], Fs);
    
    % 创建图形
    figure('Name', title_str);
    plot(f/1e6, 10*log10(Pxx), 'LineWidth', 1.2);  % 频率单位转换为MHz，功率转换为dB
    
    % 图形美化
    xlabel('频率 (MHz)', 'FontSize', 12);
    ylabel('功率谱密度 (dB/Hz)', 'FontSize', 12);
    title(title_str, 'FontSize', 14);
    grid on;
    box on;
    xlim([0 Fs/2/1e6]);  % 显示到奈奎斯特频率
    set(gca, 'FontSize', 10);
    
end
