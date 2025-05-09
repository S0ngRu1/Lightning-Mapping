%% Step1 读取引雷点的二维定位结果（需要条件筛选出合格的）
% 引入变量：位置，方位角，仰角
chj_signal_length = 5120;
match_signal_length = 6000;
yld_result_path = 'result_yld_window5120_3e8.txt';
start_signal_loc = 3.9e8;
end_signal_loc = 4e8;
step = 1e7;
% 引入两个站的位置关系
yld_sit = [0, 0, 0];
chj_sit = [1991, -7841.2, 0];
% yld相对于chj的位置
p = chj_sit-yld_sit;
dist = 8.09e3; %单位：米
c = 0.299792458;
W = 30000; % 时间误差
offsets = [-85000, -79000, -75000, -65000, -57000, -48000, -36000, -30000, -25000, -14000, -5000, 5000, 11000, 16000, 26000, 34000, 38000];
% offsets = -65000;
signal_length=1e7;
% 所有信号的开始位置
all_start_signal_loc = start_signal_loc:step:end_signal_loc;
S_results = [];
for j = 1:numel(all_start_signal_loc)-1
    start_read_loc_yld = all_start_signal_loc(j);
    end_read_loc_yld = all_start_signal_loc(j+1);
    % 记录处理的位置
    fprintf('正在处理的信号位置：%d -- %d \n', start_read_loc_yld, end_read_loc_yld);
    [yld_start_loc, yld_azimuth, yld_elevation, yld_Rcorr, yld_t123] = read_result(yld_result_path,start_read_loc_yld, end_read_loc_yld);
    yld_ch1 =read_signal('..\\20240822165932.6610CH1.dat',signal_length,start_read_loc_yld);
    chj_ch1 =read_signal('..\\2024 822 85933.651462CH1.dat',signal_length,start_read_loc_yld+ 34151156 - offsets(j));
    chj_ch2 =read_signal('..\\2024 822 85933.651462CH2.dat',signal_length,start_read_loc_yld+ 34151156 - offsets(j));
    chj_ch3 =read_signal('..\\2024 822 85933.651462CH3.dat',signal_length,start_read_loc_yld+ 34151156 - offsets(j) +165/5);
    filtered_yld_signal1 = filter_bp(yld_ch1,30e6,80e6,5);
    filtered_chj_signal1 = filter_bp(chj_ch1,30e6,80e6,5);
    filtered_chj_signal2 = filter_bp(chj_ch2,30e6,80e6,5);
    filtered_chj_signal3 = filter_bp(chj_ch3,30e6,80e6,5);
    S_results = [];
    match_results = struct('yld_start_loc', {}, 'chj_loc', {}, 'r_gccs', {}, 'dlta',{});
    h = waitbar(0, 'Processing...');
    %% Step2 根据引雷点的信号窗口得到匹配到的从化局的信号
    for i =1 :numel(yld_start_loc)
        sub_S_results = [];
        sub_R_gccs = [];
        dltas = [];
        sub_chj_locs = [];
        waitbar(i/numel(yld_start_loc), h, sprintf('位置：%d -- %d ；进度： %.2f%%', start_read_loc_yld, end_read_loc_yld, i/numel(yld_start_loc)*100));
        if yld_Rcorr(i) < 0.3 && yld_t123(i) > 1
            continue
        end
        % 转换绝对位置到相对位置
        yld_signal_start_loc = yld_start_loc(i) - start_read_loc_yld;
        start_read_loc_chj = yld_signal_start_loc;
        if yld_signal_start_loc+chj_signal_length > signal_length || yld_signal_start_loc + 1 < 0
            continue
        end
        processed_yld_signal = filtered_yld_signal1(yld_signal_start_loc+1 : yld_signal_start_loc+chj_signal_length);
        processed_yld_signal = real(windowsignal(detrend(processed_yld_signal)));

        % 读取从化局与引雷点同位置前后加6000点的信号
        if start_read_loc_chj-match_signal_length+1 <= 0 || start_read_loc_chj+match_signal_length+chj_signal_length > signal_length
            continue;
        end
        chj_match_signal1 = filtered_chj_signal1(start_read_loc_chj-match_signal_length+1:start_read_loc_chj+match_signal_length+chj_signal_length);
        chj_match_signal2 = filtered_chj_signal2(start_read_loc_chj-match_signal_length+1:start_read_loc_chj+match_signal_length+chj_signal_length);
        chj_match_signal3 = filtered_chj_signal3(start_read_loc_chj-match_signal_length+1:start_read_loc_chj+match_signal_length+chj_signal_length);

        % 寻找峰值
        [peaks, locs] = findpeaks(chj_match_signal1, 'MinPeakHeight', 13, 'MinPeakDistance', 512);
        all_locs = locs;
        % 遍历所有峰值
        num_peaks = numel(all_locs);
        if num_peaks == 0
            continue;
        end

        for pi = 1:num_peaks
            idx = all_locs(pi);
            if idx - (chj_signal_length / 2 - 1) <= 0 || idx + (chj_signal_length / 2) > match_signal_length*2+chj_signal_length
                continue;
            end
            processed_chj_signal1 = chj_match_signal1(idx - (chj_signal_length / 2)+ 1:idx + (chj_signal_length / 2));
            processed_chj_signal2 = chj_match_signal2(idx - (chj_signal_length / 2)+ 1:idx + (chj_signal_length / 2));
            processed_chj_signal3 = chj_match_signal3(idx - (chj_signal_length / 2)+ 1:idx + (chj_signal_length / 2));

            processed_chj_signal1 = real(windowsignal(detrend(processed_chj_signal1)));
            processed_chj_signal2 = real(windowsignal(detrend(processed_chj_signal2)));
            processed_chj_signal3 = real(windowsignal(detrend(processed_chj_signal3)));

            [r_gcc, lags_gcc] = xcorr(processed_yld_signal, processed_chj_signal1, 'normalized');
            R_gcc = max(r_gcc);
            t_gcc = cal_tau(r_gcc, lags_gcc');
            chj_start_idx = idx - (chj_signal_length / 2)+ 1 - t_gcc;
            chj_end_idx  = idx + (chj_signal_length / 2) - t_gcc;
            if chj_end_idx > length(chj_match_signal1) || chj_end_idx < 1 || chj_start_idx < 1 || chj_start_idx > length(chj_match_signal1)
                continue
            end
            processed_chj_signal1 = chj_match_signal1(chj_start_idx:chj_end_idx);
            processed_chj_signal2 = chj_match_signal2(chj_start_idx:chj_end_idx);
            processed_chj_signal3 = chj_match_signal3(chj_start_idx:chj_end_idx);

            processed_chj_signal1 = real(windowsignal(detrend(processed_chj_signal1)));
            processed_chj_signal2 = real(windowsignal(detrend(processed_chj_signal2)));
            processed_chj_signal3 = real(windowsignal(detrend(processed_chj_signal3)));
            [chj_start_loc, chj_azimuth, chj_elevation, chj_Rcorr, chj_t123] = get_2d_result_single_window(start_read_loc_chj,processed_chj_signal1,processed_chj_signal2,processed_chj_signal3,'chj');
            if chj_start_loc == 0 && chj_Rcorr < 0.3 && chj_t123 > 1
                continue
            end
            [R1_x, R1_y, R1_z] = sph2cart(deg2rad(90-yld_azimuth(i)), deg2rad(yld_elevation(i)),1);
            [R2_x, R2_y, R2_z] = sph2cart(deg2rad(90-chj_azimuth), deg2rad(chj_elevation),1);
%             p = chj_sit(:) - yld_sit(:);
%             A1 = [R1_x, R1_y, R1_z]';
%             A2 = [R2_x, R2_y, R2_z]';
%             % 计算叉乘并检查共线性
%             C = cross(A1, A2);
%             if norm(C) < eps
%                 continue;  % 或处理共线情况
%             end
%             c_unit = C / norm(C);
%             
%             % 构建矩阵 M = [A1, -A2, c_unit]
%             M = [A1, -A2, c_unit];
%             
%             % 使用矩阵除法替代克莱姆法则
%             solution = M \ p;
%             R1_value = solution(1);
%             R2_value = solution(2);
%             R3_value = solution(3);
%             
%             % 计算定位点（取两站解的平均值）
%             S_yld = yld_sit(:) + R1_value * A1;
%             S_chj = chj_sit(:) + R2_value * A2;
%             sub_S = (S_yld + S_chj) / 2;
%             sub_S = sub_S';

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
                sub_S = R1 + (R1_value / R2_value)*(R1_value / (R1_value + R2_value)) * R3;
            else
                % 使用第二个公式
                sub_S = R2 - (R2_value / R1_value)* (R2_value / (R1_value + R2_value))  * R3 + p;
            end
            if ~isempty(sub_S)
                t_chj = sqrt(sum((sub_S - chj_sit).^2))/c;
                t_yld = sqrt(sum((sub_S - yld_sit).^2))/c;
                dlta_t = abs(t_yld-t_chj);
                dlta_T = abs(t_gcc)*5;
                dlta = abs(dlta_t-dlta_T);
                if dlta <= W
                    dltas = [dltas;dlta];
                    sub_S_results = [sub_S_results; sub_S];
                    sub_R_gccs = [sub_R_gccs;R_gcc];
                    sub_chj_locs = [sub_chj_locs;chj_start_idx];
                end
            end
        end
        [max_R_gcc, max_R_gcc_index] = max(sub_R_gccs);
        S_results = [S_results;sub_S_results(max_R_gcc_index,:)];
        match_results = [match_results; struct('yld_start_loc', yld_start_loc(i), 'chj_loc', sub_chj_locs(max_R_gcc_index)+ start_read_loc_yld + 34151156 - offsets(j) +start_read_loc_chj-match_signal_length+1, 'r_gccs', max_R_gcc, 'dlta',dltas(max_R_gcc_index))];
    end
    close(h);
    % 创建文件名
    S_results = [S_results; S_results];
    match_results_filename = sprintf('%d_%d_match_results.mat', start_read_loc_yld, end_read_loc_yld);
    S_results_filename = sprintf('%d_%d_S_results.mat', start_read_loc_yld, end_read_loc_yld);
    % 保存变量到文件
    save(match_results_filename, 'match_results');
    save(S_results_filename, 'S_results');
end



%% Step 5: 差分到达时间 (DTOA) 技术


% 绘制 S 的结果 设置过滤条件
x_range = [-50000, 50000]; % X 的合理范围
y_range = [-50000, 50000]; % Y 的合理范围
z_range = [0, 5000];    % Z 的合理范围（Z > 0）

% 过滤数据
filtered_S = S_results(...
    S_results(:,1) >= x_range(1) & S_results(:,1) <= x_range(2) & ... % X 在合理范围内
    S_results(:,2) >= y_range(1) & S_results(:,2) <= y_range(2) & ... % Y 在合理范围内
    S_results(:,3) >= z_range(1) & S_results(:,3) <= z_range(2), :); % Z 在合理范围内


% 绘制过滤后的数据
figure;
scatter3(filtered_S(:,1), filtered_S(:,2), filtered_S(:,3),3, 'filled');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('过滤后的所有 S 点的分布');
grid on;


x = filtered_S(:, 1);
y = filtered_S(:, 2);
z = filtered_S(:, 3);

% --- 绘制二维图像 ---

num_points = size(filtered_S, 1); % 定位结果的数量

% 预分配存储角度的数组
azimuths = zeros(num_points, 1);
elevations = zeros(num_points, 1);

% 遍历每一个定位结果点
for i = 1:num_points
    % 当前定位结果点 (确保是行向量)
    current_S = filtered_S(i, :);

    % 计算从零点位置指向当前定位结果点的向量
    direction_vector = current_S(:) - yld_sit(:);

    % 调用 cart2sph_standard 函数计算方位角和仰角
    [azimuths(i), elevations(i)] = cart2sph_standard(direction_vector);
end

% 方位角-仰角散点图
figure;
scatter(azimuths, elevations, 1, 'filled');
xlabel('方位角 (度)');
ylabel('仰角 (度)');
title(['从站点 ', num2str(yld_sit), ' 看闪电的方位角 vs 仰角']);
grid on;
xlabel('Azimuth');
xlim([0, 360]);
xticks(0:40:360);
ylabel('Elevation');
ylim([-40, 100]);
yticks(-40:20:100);
