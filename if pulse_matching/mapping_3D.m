%% Step1 读取引雷点的二维定位结果（需要条件筛选出合格的）
% 引入变量：位置，方位角，仰角
chj_signal_length = 1024;
yld_result_path = 'result_yld3.5-5.5.txt';
start_read_loc_yld = 450148968;
end_read_loc_yld = 549929198;    % 引入两个站的位置关系
yld_sit = [0, 0, 0];
chj_sit = [7.8115e3, 2.1045e3, 0];
% yld相对于chj的位置
p = chj_sit-yld_sit;
dist = 8.09e3; %单位：米
c = 0.299552816;
% W = 40; % 时间误差，单位：采样点
S_results = [];
[yld_start_loc, yld_azimuth, yld_elevation, yld_Rcorr, yld_t123] = read_result(yld_result_path,start_read_loc_yld, end_read_loc_yld);
yld_azimuth = mod(yld_azimuth - 90, 360);
h = waitbar(0, 'Processing...');
%% Step2 根据引雷点的信号窗口得到匹配到的从化局的信号

for i =1 :numel(yld_start_loc)
    waitbar(i/numel(yld_start_loc), h, sprintf('Processing %.2f%%', i/numel(yld_start_loc)*100));
    if yld_Rcorr(i) < 0.3 && yld_t123(i) > 1
        continue
    end
    %     [start_read_loc_chj_top10, top10_r_gccs] =  get_match_single_yld_chj(yld_start_loc(i));
    start_read_loc_chj = yld_signal_start_loc-7727554;
    %     for j = 1:numel(start_read_loc_chj_top10)
    chj_signal1 = read_signal('../2024 822 85933.651462CH1.dat',chj_signal_length,start_read_loc_chj);
    chj_signal2 = read_signal('../2024 822 85933.651462CH2.dat',chj_signal_length,start_read_loc_chj);
    chj_signal3 = read_signal('../2024 822 85933.651462CH3.dat',chj_signal_length,start_read_loc_chj+165/5);
    [chj_start_loc, chj_azimuth, chj_elevation, chj_Rcorr, chj_t123] = get_2d_result_single_window(start_read_loc_chj,chj_signal1,chj_signal2,chj_signal3);
    if chj_start_loc == 0
        continue
    end
    chj_azimuth = mod(chj_azimuth - 90, 360);
    [R1_x, R1_y, R1_z] = az_el_to_direction(yld_azimuth(i), yld_elevation(i));
    [R2_x, R2_y, R2_z] = az_el_to_direction(chj_azimuth, chj_elevation);
    A1 = [R1_x, R1_y, R1_z];
    A2 = [R2_x, R2_y, R2_z];
    C = cross(A1, A2);
    if C == [0, 0, 0]
        continue
    end
    M = [A1; A2; C];
    % 使用克莱姆法则求R1,R2,R3的标量
    [R1_value, R2_value, R3_value] = cramer_rule(M, p);
    R1 = R1_value * A1;
    R2 = R2_value * A2;
    R3 = R3_value/norm(C)* C;
    if R1_value <= R2_value
        % 使用第一个公式
        S = R1 + (R1_value / R2_value)*(R1_value / (R1_value + R2_value)) * (R1_value / R2_value) * R3;
    else
        % 使用第二个公式
        S = R2 - (R2_value / R1_value)* (R2_value / (R1_value + R2_value)) * (R2_value / R1_value) * R3 + p;
    end
    if ~isempty(S)
        S_results = [S_results; S];
    end

    %         t_chj = sqrt(sum((S - chj_sit).^2))/c;
    %         t_yld = sqrt(sum((S - yld_sit).^2))/c;
    %         dlta_t = abs(t_yld-t_chj);
    %         dlta_T = abs(start_read_loc_chj_top10(j)-yld_start_loc(i))/5;
    %         if abs(dlta_t-dlta_T) <= W
    %         end
    %% Step 4: 保存S 并计算源到两个站点的时间延迟
end
close(h);

%% Step 5: 差分到达时间 (DTOA) 技术


% 绘制 S 的结果
% 设置过滤条件
x_range = [-500000, 500000]; % X 的合理范围
y_range = [-500000, 500000]; % Y 的合理范围
z_range = [0, 500000];    % Z 的合理范围（Z > 0）

% 过滤数据
filtered_S = S_results(...
    S_results(:,1) >= x_range(1) & S_results(:,1) <= x_range(2) & ... % X 在合理范围内
    S_results(:,2) >= y_range(1) & S_results(:,2) <= y_range(2) & ... % Y 在合理范围内
    S_results(:,3) >= z_range(1) & S_results(:,3) <= z_range(2), :); % Z 在合理范围内


% 绘制过滤后的数据
figure;
scatter3(filtered_S(:,1), filtered_S(:,2), filtered_S(:,3), 1, 'filled');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('过滤后的所有 S 点的分布');
grid on;

% 转换为极坐标
[x, y, z] = deal(filtered_S(:,1), filtered_S(:,2), filtered_S(:,3));
theta = atan2(y, x); % 计算角度 (弧度)
r = sqrt(x.^2 + y.^2); % 计算半径

% 绘制极坐标图
figure;
polarscatter(theta, r, 1, z, 'filled'); % 使用颜色表示 Z 坐标
title('过滤后的所有 S 点的极坐标分布');
colorbar; % 添加颜色条表示 Z 值