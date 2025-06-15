% --- 1. 定义已知参数 ---
yld_sit = [0, 0, 0];         % 观测站1 (B1) 的位置 (x,y,z)
chj_sit = [1991, -7841.2, -27]; % 观测站2 (B2) 的位置 (x,y,z)

% 定义一个已知的目标三维点
target_true_pos = [15000, -30000, 2000]; % 例如 (x,y,z)

fprintf('--- 输入参数 ---\n');
fprintf('观测站1 (YLD) 位置: [%.2f, %.2f, %.2f]\n', yld_sit(1), yld_sit(2), yld_sit(3));
fprintf('观测站2 (CHJ) 位置: [%.2f, %.2f, %.2f]\n', chj_sit(1), chj_sit(2), chj_sit(3));
fprintf('真实目标位置:       [%.2f, %.2f, %.2f]\n\n', target_true_pos(1), target_true_pos(2), target_true_pos(3));


[yld_azimuth_true_deg, yld_elevation_true_deg] = calculate_angles_for_algorithm(yld_sit, target_true_pos);
[chj_azimuth_true_deg, chj_elevation_true_deg] = calculate_angles_for_algorithm(chj_sit, target_true_pos);

fprintf('--- 计算得到的输入角度 ---\n');
fprintf('YLD站到目标的方位角 (算法输入): %.4f deg\n', yld_azimuth_true_deg);
fprintf('YLD站到目标的俯仰角 (算法输入): %.4f deg\n', yld_elevation_true_deg);
fprintf('CHJ站到目标的方位角 (算法输入): %.4f deg\n', chj_azimuth_true_deg);
fprintf('CHJ站到目标的俯仰角 (算法输入): %.4f deg\n\n', chj_elevation_true_deg);


% 调用算法
estimated_pos = localization_algorithm(yld_sit, chj_sit, ...
                                       yld_azimuth_true_deg, yld_elevation_true_deg, ...
                                       chj_azimuth_true_deg, chj_elevation_true_deg);

% --- 4. 结果对比 ---
fprintf('--- 定位结果 ---\n');
fprintf('算法估算的目标位置: [%.4f, %.4f, %.4f]\n', estimated_pos(1), estimated_pos(2), estimated_pos(3));

if all(~isnan(estimated_pos))
    error_vec = estimated_pos - target_true_pos;
    error_magnitude = norm(error_vec);
    fprintf('定位误差向量:       [%.2e, %.2e, %.2e]\n', error_vec(1), error_vec(2), error_vec(3));
    fprintf('定位误差大小:       %.2e\n', error_magnitude);

    if error_magnitude < 1e-9 % 根据期望的精度调整阈值
        fprintf('\n结论: 算法结果与真实目标位置非常接近，验证通过！\n');
    else
        fprintf('\n结论: 算法结果与真实目标位置存在差异，请检查。\n');
    end
else
    fprintf('\n结论: 算法未能成功计算目标位置。\n');
end



% --- 2. 正向计算：从站点到目标的真实方位角和俯仰角 ---


% 函数：根据站点和目标位置，计算算法所需的方位角和俯仰角
% 注意：此处的方位角是您算法 sph2cart(deg2rad(90-THIS_AZIMUTH), ...) 中的 THIS_AZIMUTH
function [algo_az_deg, algo_el_deg] = calculate_angles_for_algorithm(station_pos, target_pos)
    delta_vec = target_pos - station_pos;
    dx = delta_vec(1);
    dy = delta_vec(2);
    dz = delta_vec(3);

    % cart2sph 返回的方位角是从x轴到点在xy平面投影的逆时针角度
    % 返回的俯仰角是从xy平面到点的角度
    [az_cart_rad, el_cart_rad, ~] = cart2sph(dx, dy, dz);

    az_cart_deg = rad2deg(az_cart_rad); % 标准的从x轴计算的方位角 (°)
    el_cart_deg = rad2deg(el_cart_rad); % 标准的从xy平面计算的俯仰角 (°)

    % 算法内部使用 sph2cart(deg2rad(90-input_az), deg2rad(input_el), 1)
    % 这意味着 (90-input_az) 对应于 az_cart_deg
    % 所以，input_az = 90 - az_cart_deg
    algo_az_deg = 90.0 - az_cart_deg; % 使用浮点数确保精度
    algo_el_deg = el_cart_deg;
end




% --- 3. 定位算法 ---
function estimated_target_global_pos = localization_algorithm(station1_pos, station2_pos, s1_az_deg, s1_el_deg, s2_az_deg, s2_el_deg)
    p_baseline_vec = station2_pos - station1_pos;
    % 从球坐标获取单位方向向量 A1, A2 
    [R1_x, R1_y, R1_z] = sph2cart(deg2rad(90-s1_az_deg), deg2rad(s1_el_deg),1);
    A1 = [R1_x, R1_y, R1_z];

    [R2_x, R2_y, R2_z] = sph2cart(deg2rad(90-s2_az_deg), deg2rad(s2_el_deg),1);
    A2 = [R2_x, R2_y, R2_z];

    C = cross(A1, A2);
    norm_C = norm(C);

    if norm_C < eps % 使用 MATLAB 的 eps 来处理浮点数比较
        warning('方向向量 A1 和 A2 几乎平行，无法准确定位。');
        estimated_target_global_pos = [NaN, NaN, NaN];
        return;
    end
    c_unit = C / norm_C;

    M = [A1(1), -A2(1), c_unit(1);
         A1(2), -A2(2), c_unit(2);
         A1(3), -A2(3), c_unit(3)];

    detM = det(M);
    if abs(detM) < eps % 检查行列式是否过小
        warning('矩阵 M 接近奇异或病态，结果可能不可靠。');
        estimated_target_global_pos = [NaN, NaN, NaN];
        return;
    end

    % 使用克莱姆法则求解 R1_value, R2_value, R3_value (标量)
    % p_baseline_vec 需要是列向量以便与M的列进行替换
    p_col_vec = p_baseline_vec';

    detR1 = det([p_col_vec, M(:,2), M(:,3)]);
    detR2 = det([M(:,1), p_col_vec, M(:,3)]);
    detR3 = det([M(:,1), M(:,2), p_col_vec]);

    R1_value = detR1 / detM;
    R2_value = detR2 / detM;
    R3_value = detR3 / detM; % 这是沿 c_unit 方向的长度
    R1_vec = R1_value * A1;
    R2_vec_from_s2 = R2_value * A2; % 这是从station2出发的向量
    R3_vec = R3_value * c_unit;   % 这是S1S2向量 (R3_value/norm_C * C 等于 R3_value * c_unit)
    if R1_value <= R2_value
        factor = (R1_value / R2_value) * (R1_value / (R1_value + R2_value));
        if abs(R1_value + R2_value) < eps || abs(R2_value) < eps
             warning('R1_value + R2_value 或 R2_value 接近零，插值因子可能无效。');
             sub_S = R1_vec; 
        else
            sub_S = R1_vec + factor * R3_vec;
        end
    else
        factor = (R2_value / R1_value) * (R2_value / (R1_value + R2_value));
        if abs(R1_value + R2_value) < eps || abs(R1_value) < eps
            warning('R1_value + R2_value 或 R1_value 接近零，插值因子可能无效。');
            sub_S = p_baseline_vec + R2_vec_from_s2; 
        else
             
             sub_S = R2_vec_from_s2 - factor * R3_vec + p_baseline_vec;
        end
    end

    estimated_target_global_pos = sub_S + station1_pos;
end

