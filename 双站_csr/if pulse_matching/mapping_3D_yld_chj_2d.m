%% Step1 读取引雷点的二维定位结果（需要条件筛选出合格的）
% 引入变量：位置，方位角，仰角
yld_result_path = 'result_yld_window5120_3e8.txt';
chj_result_path = 'result_chj_window5120_newloc_all.txt';
start_signal_loc = 4.696e8;
end_signal_loc = 4.702e8;
% 引入两个站的位置关系
yld_sit = [0, 0, 0];
chj_sit = [2003.7972, -7844.7836, -27];
% yld相对于chj的位置
p = chj_sit-yld_sit;
dist = 8.0967e+03; %单位：米
c = 0.299792458;
W = 30000; % 时间误差
offsets_init = -85000;
signal_length=5120;
% 所有信号的开始位置
all_S_results = [];
all_match_results = [];
% 记录处理的位置
[yld_start_locs, yld_azimuths, yld_elevations, yld_Rcorrs, yld_t123s] = read_result(yld_result_path,start_signal_loc, end_signal_loc);
total_iterations = numel(yld_start_locs);
h_waitbar = waitbar(0, '正在处理 YLD 数据，请稍候...', 'Name', '处理进度');
for j = 1:total_iterations
    progress = j / total_iterations;
    % 更新进度条的值和文本
    waitbar(progress, h_waitbar, sprintf('正在处理 YLD 数据: %d / %d (%.1f%%)', j, total_iterations, progress * 100));

    if yld_Rcorrs(j) < 0.3 || yld_t123s(j) > 1
        continue;
    end
    yld_start_loc = yld_start_locs(j);
    chj_start_loc = yld_start_loc + 34151156 - offsets_init - floor((yld_start_loc-start_signal_loc)/1282) - signal_length - 6000;
    chj_end_loc = chj_start_loc + signal_length*2 + 2*6000;
    [chj_start_locs, chj_azimuths, chj_elevations, chj_Rcorrs, chj_t123s] = read_result(chj_result_path,chj_start_loc, chj_end_loc);
    yld_ch1 =read_signal('..\\20240822165932.6610CH1.dat',signal_length,yld_start_loc);
    filtered_yld_signal1 = filter_bp(yld_ch1,30e6,80e6,5);
    sub_Rcorrs = [];
    for k = 1:numel(chj_start_locs)
        if chj_Rcorrs(k) <  0.2 || chj_t123s(k) > 1
            sub_Rcorrs =[sub_Rcorrs; 0];
            continue;
        end
        chj_ch1 =read_signal('..\\2024 822 85933.651462CH1.dat', signal_length, chj_start_locs(k));
        filtered_chj_signal1 = filter_bp(chj_ch1,30e6,80e6,5);
        [r_gcc, lags_gcc] = xcorr(filtered_yld_signal1, filtered_chj_signal1, 'normalized');
        sub_Rcorrs =[sub_Rcorrs; max(r_gcc)];
    end
    [max_R_gcc, max_R_gcc_index] = max(sub_Rcorrs);
    if max_R_gcc < 0.08
        continue;
    end
    chj_start_loc_max = chj_start_locs(max_R_gcc_index);
    chj_azimuth = chj_azimuths(max_R_gcc_index);
    chj_elevation = chj_elevations(max_R_gcc_index);
    S_result = [];
    match_result = struct('yld_start_loc', {}, 'chj_start_loc', {}, 'r_gcc', {}, 'chj_azimuth',{},'chj_elevation',{},'dlta',{});
    [R1_x, R1_y, R1_z] = sph2cart(deg2rad(90-yld_azimuths(j)), deg2rad(yld_elevations(j)),1);
    [R2_x, R2_y, R2_z] = sph2cart(deg2rad(90-chj_azimuth), deg2rad(chj_elevation),1);
    A1 = [R1_x, R1_y, R1_z];
    A2 = [R2_x, R2_y, R2_z];
    C = cross(A1, A2);
    if norm(c) < eps
        continue;  % 避免除以零
    end
    c_unit = C  / norm(C);  % 单位向量
    M = [A1(1), -A2(1), c_unit(1);
        A1(2), -A2(2), c_unit(2);
        A1(3), -A2(3), c_unit(3)];
    % 使用克莱姆法则求R1,R2,R3的标量
    detM = det(M);
    detR1 = det([p', M(:,2), M(:,3)]);
    detR2 = det([M(:,1), p', M(:,3)]);
    detR3 = det([M(:,1), M(:,2), p']);
    R1_value = detR1 / detM;
    R2_value = detR2 / detM;
    R3_value = detR3 / detM;
    R1 = R1_value * A1;
    R2 = R2_value * A2;
    R3 = R3_value/norm(C)* C;
    if R1_value <= R2_value
        % 使用第一个公式
        S = R1 + (R1_value / R2_value)*(R1_value / (R1_value + R2_value)) * R3;
    else
        % 使用第二个公式
        S = R2 - (R2_value / R1_value)* (R2_value / (R1_value + R2_value))  * R3 + p;
    end
    if S(3) < 0
        S = -S;
    end
    if ~isempty(S)
        t_chj = sqrt(sum((S - chj_sit).^2))/c;
        t_yld = sqrt(sum((S - yld_sit).^2))/c;
        dlta_t = abs(t_yld-t_chj);
        dlta_T = abs(chj_start_loc_max-chj_start_loc- signal_length - 6000)*5;
        dlta = abs(dlta_t-dlta_T);
        if dlta <= W
            all_S_results = [all_S_results; S];
            all_match_results = [all_match_results; struct('yld_start_loc', yld_start_locs(j), 'chj_loc', chj_start_loc_max,'chj_azimuth',chj_azimuth,'chj_elevation',chj_elevation, 'r_gccs', max_R_gcc, 'dlta',dlta)];
        end
    end
end
close(h_waitbar);
plot_3d;
