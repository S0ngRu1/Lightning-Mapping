clear; clc; 
% close all;
yld_result_path = 'result_yld_4.6-5.1e8_th(15)_bp.txt';
start_read_loc_yld = 469400180;
end_read_loc_yld = 470019709;  
signal_length = 1e7;
r_loction = 4.693e8;
yld_ch1 =read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction);
chj_ch1 =read_signal('..\\2024 822 85933.651462CH1.dat',signal_length,r_loction+ 34187200);
filtered_yld_signal1 = filter_bp(yld_ch1,30e6,80e6,5);
filtered_chj_signal1 = filter_bp(chj_ch1,30e6,80e6,5);

window_lengths = 10240;
[yld_start_loc, yld_azimuth, yld_elevation, yld_Rcorr, yld_t123] = read_result(yld_result_path,start_read_loc_yld, end_read_loc_yld);
h = waitbar(0, 'Processing...');

% 初始化空结构数组
results = struct(...
    'startLoc', {}, ...           % 原始 yld 起始位置
    'matchPositions', {}, ...     % 各窗口的 CHJ 起始位置向量
    'maxCorrelations', {}, ...    % 各窗口的最大相关系数向量
    'timeDelays', {} ...          % 各窗口的时间偏差向量
);
nYlds = numel(yld_start_loc);
for j = 1:nYlds
    % 条件过滤：不满足条件的直接跳过，不加入 results
    if yld_Rcorr(j) < 0.3 && yld_t123(j) > 1
        continue;
    end
    % 更新等待条
    waitbar(j/nYlds, h, sprintf('Processing %.2f%%', j/nYlds*100));

    % 转换绝对位置到相对位置
    yld_signal_start_loc = yld_start_loc(j) - r_loction;
    current_chj_read_loc = yld_signal_start_loc;

    loc_diff = 40000;

    % 提取 yld 段并校验边界
    idxYldStart = yld_signal_start_loc - window_lengths/2 + 1;
    idxYldEnd   = yld_signal_start_loc + window_lengths/2;

    if idxYldEnd > length(filtered_yld_signal1)
        break;
    end
    yld_seg = filtered_yld_signal1(idxYldStart:idxYldEnd);
    yld_seg = real(windowsignal(detrend(yld_seg)));

    current_chj_read_loc = current_chj_read_loc - loc_diff/2;
    chj_len = window_lengths + loc_diff;
    idxChjStart = max(1,floor(current_chj_read_loc) + 1);
    idxChjEnd   = min(length(filtered_chj_signal1),floor(current_chj_read_loc + chj_len));
    if idxChjStart < 1 || idxChjEnd > length(filtered_chj_signal1) || idxChjEnd < 1
        continue;
    end
    chj_seg = filtered_chj_signal1(idxChjStart:idxChjEnd);
    chj_seg = real(windowsignal(detrend(chj_seg)));

    % 计算 GCC-PHAT 并获取峰值信息
    [r_gcc, lags_gcc] = xcorr(yld_seg, chj_seg,'none');
    t_gcc            = cal_tau(r_gcc, lags_gcc');
    r_gcc = r_gcc ./ (norm(chj_seg) * norm(yld_seg));
    [maxC, idxMax]   = max(r_gcc);
    % 记录结果
    matchPos = floor(idxChjStart);
    maxCorrs = maxC;
    timeDs   = t_gcc;

     results(end+1) = struct(...
        'startLoc',        yld_start_loc(j)- window_lengths/2 + 1, ...
        'matchPositions',  r_loction+ 34187200 + matchPos + 1, ...
        'maxCorrelations', maxCorrs, ...
        'timeDelays',      timeDs ...
    );
end



yld_ch1 =read_signal('..\\20240822165932.6610CH1.dat',10240,469400180);
yld_ch1 = filter_bp(yld_ch1,30e6,80e6,5);

chj_ch1 = read_signal('..\\2024 822 85933.651462CH1.dat',10240,503587381+(-310));
chj_ch1 = filter_bp(chj_ch1,30e6,80e6,5);
figure
subplot(2,1,1);plot(yld_ch1);title('yld');
subplot(2,1,2);plot(chj_ch1);title('chj');


chj_ch1_window =read_signal('..\\2024 822 85933.651462CH1.dat',10240+loc_diff,469415189+34187200-loc_diff/2);
chj_ch1_window = filter_bp(chj_ch1_window,30e6,80e6,5);
chj_ch1 = read_signal('..\\2024 822 85933.651462CH1.dat',10240,503587201+26208);
chj_ch1 = filter_bp(chj_ch1,30e6,80e6,5);

plot(yld_ch1,'b','LineWidth', 1.5)
hold on
plot(chj_ch1_window-200,'r','LineWidth', 1.5)

% figure
% plot(lags_gcc,r_gcc)
ch1 = read_signal('..\\20240822165932.6610CH1.dat',10240,469400180);
ch2 = read_signal('..\\20240822165932.6610CH1.dat',10240,469400180+6000);
[r_gcc, lags_gcc] = xcorr(ch1, ch2,'normalized');
t_gcc            = cal_tau(r_gcc, lags_gcc');
