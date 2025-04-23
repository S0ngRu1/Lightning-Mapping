clear; clc; close all;
yld_result_path = 'result_yld_th1_4_5e8-5_5e8.txt';
start_read_loc_yld = 469400180;
end_read_loc_yld = 470019709;  
signal_length = 8e5;
r_loction = 4.694e8;
diff = 34151156;
yld_ch1 =read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction);
chj_ch1 =read_signal('..\\2024 822 85933.651462CH1.dat',signal_length,r_loction+ diff);
win_length = 10240;
filtered_yld_signal1 = filter_bp(yld_ch1,20e6,80e6,5);
filtered_chj_signal1 = filter_bp(chj_ch1,20e6,80e6,5);
[yld_start_loc, yld_azimuth, yld_elevation, yld_Rcorr, yld_t123] = read_result(yld_result_path,start_read_loc_yld, end_read_loc_yld);
h = waitbar(0, 'Processing...');
% 批量记录每个 yld_start_loc 的匹配结果（动态添加，无预分配）
nYlds = numel(yld_start_loc);
nWins = 1;

% 初始化空结构数组
results = struct(...
    'startLoc', {}, ...           % 原始 yld 起始位置
    'matchPositions', {}, ...     % 各窗口的 CHJ 起始位置向量
    'maxCorrelations', {}, ...    % 各窗口的最大相关系数向量
    'timeDelays', {} ...          % 各窗口的时间偏差向量
);

for j = 1:nYlds
    % 条件过滤：不满足条件的直接跳过，不加入 results
    if yld_Rcorr(j) < 0.3 && yld_t123(j) > 1
        continue;
    end
    % 更新等待条
    waitbar(j/nYlds, h, sprintf('Processing %.2f%%', j/nYlds*100));

    % 转换绝对位置到相对位置
    yld_signal_start_loc = yld_start_loc(j) - 4.694e8;
    current_chj_read_loc = yld_signal_start_loc;

    % 准备记录临时向量
    matchPos    = nan(1, nWins);
    maxCorrs    = nan(1, nWins);
    timeDs      = nan(1, nWins);

    % 逐窗匹配
    for i = 1:nWins
        % 提取 yld 段并校验边界
        idxYldStart = yld_signal_start_loc + 1;
        idxYldEnd   = yld_signal_start_loc + win_length;
        if idxYldEnd > length(filtered_yld_signal1)
            break;
        end
        yld_seg = filtered_yld_signal1(idxYldStart:idxYldEnd);
        yld_seg = real(windowsignal(detrend(yld_seg)));
        % 更新 chj 读取起点并边界检查
        idxChjStart = max(1,floor(current_chj_read_loc) + 1);
        idxChjEnd   = min(length(filtered_chj_signal1),floor(current_chj_read_loc + win_length));
        if idxChjStart < 1 || idxChjEnd > length(filtered_chj_signal1) || idxChjEnd < 1
            continue;
        end
        chj_seg = filtered_chj_signal1(idxChjStart:idxChjEnd);

        % 计算 GCC-PHAT 并获取峰值信息
        [r_gcc, lags_gcc] = xcorr(chj_seg, yld_seg, 'none');
        t_gcc            = cal_tau(r_gcc, lags_gcc');
        r_gcc = r_gcc ./ (norm(chj_seg) * norm(yld_seg));
        [maxC, idxMax]   = max(r_gcc);
        % 记录本轮结果
        matchPos(i) = floor(idxChjStart);
        maxCorrs(i) = maxC;
        timeDs(i)   = t_gcc;
    end

    % 将本次计算结果追加到 results
    results(end+1) = struct(...
        'startLoc',        yld_start_loc(j), ...
        'matchPositions',  r_loction+ diff + matchPos, ...
        'maxCorrelations', maxCorrs, ...
        'timeDelays',      timeDs ...
    );
end
close(h);


yld_ch1 =read_signal('..\\20240822165932.6610CH1.dat',1e5,469467858-5e4);
chj_ch1 =read_signal('..\\2024 822 85933.651462CH1.dat',1e5,503619771+27-5e4);
subplot(2,1,1);plot(yld_ch1);title('yld');
subplot(2,1,2);plot(chj_ch1);title('chj');

