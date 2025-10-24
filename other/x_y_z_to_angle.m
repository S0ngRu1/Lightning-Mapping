A = [0, 0];
B = [31.4263,-10.0785];
C = [19.2537,-36.3713];
% 计算AB, AC, BC的距离
AB = sqrt((B(1) - A(1))^2 + (B(2) - A(2))^2);
AC = sqrt((C(1) - A(1))^2 + (C(2) - A(2))^2);
BC = sqrt((C(1) - B(1))^2 + (C(2) - B(2))^2);

% 输出三边长度
fprintf('AB = %.4f\n', AB);
fprintf('AC = %.4f\n', AC);
fprintf('BC = %.4f\n', BC);

% 计算AB, AC, BC的方向向量
v_AB = B - A; % 向量AB
v_AC = C - A; % 向量AC
v_BC = C - B; % 向量BC

% 定义正北方向的单位向量
v_north = [0, 1]; 

% 计算边AB与正北方向的夹角
cos_AB_north = dot(v_AB, v_north) / (norm(v_AB) * norm(v_north));
angle_AB_north = acos(cos_AB_north); % 返回弧度
angle_AB_north_deg = rad2deg(angle_AB_north); % 转换为度

% 计算边AC与正北方向的夹角
cos_AC_north = dot(v_AC, v_north) / (norm(v_AC) * norm(v_north));
angle_AC_north = acos(cos_AC_north); % 返回弧度
angle_AC_north_deg = rad2deg(angle_AC_north); % 转换为度

% 计算边BC与正北方向的夹角
cos_BC_north = dot(v_BC, v_north) / (norm(v_BC) * norm(v_north));
angle_BC_north = acos(cos_BC_north); % 返回弧度
angle_BC_north_deg = rad2deg(angle_BC_north); % 转换为度

% 输出结果
fprintf('Angle between AB and North: %.4f degrees\n', angle_AB_north_deg);
fprintf('Angle between AC and North: %.4f degrees\n', angle_AC_north_deg);
fprintf('Angle between BC and North: %.4f degrees\n', angle_BC_north_deg);