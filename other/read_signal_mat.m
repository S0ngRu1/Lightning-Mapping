%% 从 signal.mat 文件中读取指定片段的信号
function signal = read_signal_mat(signal_mat_path, r_length, r_location)
    % 功能：从 signal.mat 文件中读取从 r_location 开始、长度为 r_length 的信号片段
    % 输入：
    %   signal_mat_path - signal.mat 文件的完整路径（字符串）
    %   r_length        - 需要读取的信号长度（正整数，单位：采样点）
    %   r_location      - 起始位置（正整数，从1开始的索引，即第几个采样点）
    % 输出：
    %   signal          - 截取的信号片段（double类型，行向量）
    
    
    % 加载 signal.mat 文件
    try
        full_signal = load(signal_mat_path);  % 加载MAT文件，返回结构体
    catch ME
        error('加载文件失败：%s（请检查路径是否正确）', ME.message);
    end
    full_signal = full_signal.double_400';
    % 截取信号片段并转换为double类型（与原函数输出格式一致）
    signal = double(full_signal(r_location:r_location + r_length - 1));
end