%% GPT给的三维定位方法

P_true = [5000; 2000; 3000];  % 真实点 [m]
yld_sit = [0; 0; 0];          
chj_sit = [1991; -7841.2; 0];

% 计算站点到目标的向量
v_yld = P_true - yld_sit;
v_chj = P_true - chj_sit;    

% 计算方位角和仰角（弧度制）
Az_yld = atan2(v_yld(1), v_yld(2));                
El_yld = atan2(v_yld(3), sqrt(v_yld(1)^2 + v_yld(2)^2));  
Az_chj = atan2(v_chj(1), v_chj(2));                
El_chj = atan2(v_chj(3), sqrt(v_chj(1)^2 + v_chj(2)^2));  

% 构造方向向量
dir_yld = [sin(Az_yld)*cos(El_yld); cos(Az_yld)*cos(El_yld); sin(El_yld)];
dir_chj = [sin(Az_chj)*cos(El_chj); cos(Az_chj)*cos(El_chj); sin(El_chj)];

% 解两条射线交点: yld + t*dir_yld = chj + s*dir_chj
M = [dir_yld, -dir_chj]; 
ts = M \ (chj_sit - yld_sit);   % 求解 t 和 s
t = ts(1); 
s = ts(2);

% 计算估计点（理论上 P_est1 和 P_est2 相等）
P_est1 = yld_sit + t * dir_yld;
P_est2 = chj_sit + s * dir_chj;
P_est = 0.5*(P_est1 + P_est2);  % 平均化以降低数值误差
fprintf("----------------------GPT给的三维定位方法----------------------\n")
fprintf('真实位置: [%.2f, %.2f, %.2f]\n', P_true);
fprintf('估计位置: [%.2f, %.2f, %.2f]\n', P_est);
fprintf('定位误差(米): %.6f\n', norm(P_est - P_true));


clear
%% 咱们的三维定位方法
fprintf("----------------------咱们的三维定位方法----------------------\n")
s_true = [5000; 2000; 3000];  % 单位米
yld_sit = [0, 0, 0];
chj_sit = [2003.7972, -7844.7836, -27];
p = chj_sit - yld_sit;
% 站点到真实点的向量
v_yld = s_true - yld_sit;
v_chj = s_true - chj_sit;

% 站点yld的方位角、仰角
Az_yld = atan2(v_yld(1), v_yld(2));                % 方位角：东为0°，逆时针
El_yld = atan2(v_yld(3), sqrt(v_yld(1)^2 + v_yld(2)^2));  % 仰角
% Az_yld = atan2(v_yld(2),-v_yld(1));                 % 注意x取负值，因为x轴正方向为南
% Az_yld = rad2deg(Az_yld);
% El_yld = asin(v_yld(3)/sqrt(v_yld(1)^2 + v_yld(2)^2+ v_yld(3)^2));  % 仰角
% El_yld = rad2deg(El_yld);

% 站点chj的方位角、仰角
Az_chj = atan2(v_chj(1), v_chj(2));                
El_chj = atan2(v_chj(3), sqrt(v_chj(1)^2 + v_chj(2)^2));  
% Az_chj = atan2(v_chj(2),-v_chj(1));                 % 注意x取负值，因为x轴正方向为南
% Az_chj = rad2deg(Az_chj);
% El_chj = asin(v_chj(3)/sqrt(v_chj(1)^2 + v_chj(2)^2+ v_chj(3)^2));  % 仰角
% El_chj = rad2deg(El_chj);
% 转成角度
yld_azimuth = rad2deg(Az_yld);
yld_elevation = rad2deg(El_yld);
chj_azimuth = rad2deg(Az_chj);
chj_elevation = rad2deg(El_chj);


% 构造方向向量
[R1_x, R1_y, R1_z] = sph2cart(deg2rad(90 - yld_azimuth), deg2rad(yld_elevation), 1);
[R2_x, R2_y, R2_z] = sph2cart(deg2rad(90 - chj_azimuth), deg2rad(chj_elevation), 1);

A1 = [R1_x, R1_y, R1_z];
A2 = [R2_x, R2_y, R2_z];

% 叉乘求公共法向量
C = cross(A1, A2);
if norm(C) < eps
    error('方向向量共线，无法定位');
end
c_unit = C / norm(C);  % 单位向量

% 构建矩阵
M = [A1(1), -A2(1), c_unit(1);
     A1(2), -A2(2), c_unit(2);
     A1(3), -A2(3), c_unit(3)];

% 克莱姆法则
detM = det(M);
detR1 = det([p', M(:,2), M(:,3)]);
detR2 = det([M(:,1), p', M(:,3)]);
detR3 = det([M(:,1), M(:,2), p']);

R1_value = detR1 / detM;
R2_value = detR2 / detM;
R3_value = detR3 / detM;

R1 = R1_value * A1;
R2 = R2_value * A2;
R3 = R3_value / norm(C) * C;
R3 = -R3;  
% 计算交点 sub_S
if R1_value <= R2_value
    sub_S = R1 + (R1_value / R2_value)*(R1_value / (R1_value + R2_value)) * R3;
else
    sub_S = R2 - (R2_value / R1_value)*(R2_value / (R1_value + R2_value)) * R3 + p;
end
% 还原到地面坐标
S_result = sub_S + yld_sit;  % 这里注意，需要加上yld_sit作为起点偏移
% 打印真实点 vs 估计点
fprintf('真实位置: [%.2f, %.2f, %.2f]\n', s_true);
fprintf('估计位置: [%.2f, %.2f, %.2f]\n', S_result);
fprintf('定位误差(米): %.6f\n', norm(S_result - s_true));



clear
%% 修改方位角转方向向量的三维定位方法
fprintf("----------------------修改方位角转方向向量的三维定位方法----------------------\n")
P_true = [5000; 2000; 3000];  % 单位米
yld_sit = [0, 0, 0];
chj_sit = [1991, -7841.2, 0];
p = chj_sit - yld_sit;
% 站点到真实点的向量
v_yld = P_true - yld_sit;
v_chj = P_true - chj_sit;

% 站点yld的方位角、仰角
Az_yld = atan2(v_yld(1), v_yld(2));                % 方位角：东为0°，逆时针
El_yld = atan2(v_yld(3), sqrt(v_yld(1)^2 + v_yld(2)^2));  % 仰角

% 站点chj的方位角、仰角
Az_chj = atan2(v_chj(1), v_chj(2));                
El_chj = atan2(v_chj(3), sqrt(v_chj(1)^2 + v_chj(2)^2));  

% 构造方向向量
A1 = [sin(Az_yld)*cos(El_yld); cos(Az_yld)*cos(El_yld); sin(El_yld)]';
A2 = [sin(Az_chj)*cos(El_chj); cos(Az_chj)*cos(El_chj); sin(El_chj)]';

% 叉乘求公共法向量
C = cross(A1, A2);
if norm(C) < eps
    error('方向向量共线，无法定位');
end
c_unit = C / norm(C);  % 单位向量

% 构建矩阵
M = [A1(1), -A2(1), c_unit(1);
     A1(2), -A2(2), c_unit(2);
     A1(3), -A2(3), c_unit(3)];

% 克莱姆法则
detM = det(M);
detR1 = det([p', M(:,2), M(:,3)]);
detR2 = det([M(:,1), p', M(:,3)]);
detR3 = det([M(:,1), M(:,2), p']);

R1_value = detR1 / detM;
R2_value = detR2 / detM;
R3_value = detR3 / detM;

R1 = R1_value * A1;
R2 = R2_value * A2;
R3 = R3_value / norm(C) * C;
R3 = -R3;  
% 计算交点 sub_S
if R1_value <= R2_value
    sub_S = R1 + (R1_value / R2_value)*(R1_value / (R1_value + R2_value)) * R3;
else
    sub_S = R2 - (R2_value / R1_value)*(R2_value / (R1_value + R2_value)) * R3 + p;
end
% 还原到地面坐标
S_result = sub_S + yld_sit;  % 这里注意，需要加上yld_sit作为起点偏移
% 打印真实点 vs 估计点
fprintf('真实位置: [%.2f, %.2f, %.2f]\n', P_true);
fprintf('估计位置: [%.2f, %.2f, %.2f]\n', S_result);
fprintf('定位误差(米): %.6f\n', norm(S_result - P_true));


%% 修改三维定位公式的方法
fprintf("----------------------修改三维定位公式的方法----------------------\n")
P_true = [5000; 2000; 3000];  % 单位米
yld_sit = [0, 0, 0]';
chj_sit = [1991, -7841.2, 0]';
p = chj_sit - yld_sit;
% 站点到真实点的向量
v_yld = P_true - yld_sit;
v_chj = P_true - chj_sit;

% 站点yld的方位角、仰角
Az_yld = atan2(v_yld(1), v_yld(2));                % 方位角：东为0°，逆时针
El_yld = atan2(v_yld(3), sqrt(v_yld(1)^2 + v_yld(2)^2));  % 仰角

% 站点chj的方位角、仰角
Az_chj = atan2(v_chj(1), v_chj(2));                
El_chj = atan2(v_chj(3), sqrt(v_chj(1)^2 + v_chj(2)^2));  

% 转成角度
yld_azimuth = rad2deg(Az_yld);
yld_elevation = rad2deg(El_yld);
chj_azimuth = rad2deg(Az_chj);
chj_elevation = rad2deg(El_chj);


% 构造方向向量
[R1_x, R1_y, R1_z] = sph2cart(deg2rad(90 - yld_azimuth), deg2rad(yld_elevation), 1);
[R2_x, R2_y, R2_z] = sph2cart(deg2rad(90 - chj_azimuth), deg2rad(chj_elevation), 1);

dir_yld = [R1_x, R1_y, R1_z]';
dir_chj = [R2_x, R2_y, R2_z]';


% 解两条射线交点: yld + t*dir_yld = chj + s*dir_chj
M = [dir_yld, -dir_chj]; 
ts = M \ (chj_sit - yld_sit);   % 求解 t 和 s
t = ts(1); 
s = ts(2);

% 计算估计点（理论上 P_est1 和 P_est2 相等）
P_est1 = yld_sit + t * dir_yld;
P_est2 = chj_sit + s * dir_chj;
P_est = 0.5*(P_est1 + P_est2);  % 平均化以降低数值误差
fprintf('真实位置: [%.2f, %.2f, %.2f]\n', P_true);
fprintf('估计位置: [%.2f, %.2f, %.2f]\n', P_est);
fprintf('定位误差(米): %.6f\n', norm(P_est - P_true));



%%
%%验证转方向向量
fprintf("----------------------验证方位角仰角转方向向量的方法----------------------\n")
P_true = [5000, 2000, 3000];
yld_sit = [0, 0, 0];
chj_sit = [1991, -7841.2, 0];
% 计算站点到目标的向量
v_yld = P_true - yld_sit;
v_chj = P_true - chj_sit; 

% 计算方位角和仰角（弧度制）
Az_yld = atan2(v_yld(1), v_yld(2));                
El_yld = atan2(v_yld(3), sqrt(v_yld(1)^2 + v_yld(2)^2));  
Az_chj = atan2(v_chj(1), v_chj(2));                
El_chj = atan2(v_chj(3), sqrt(v_chj(1)^2 + v_chj(2)^2));  

% 构造方向向量
dir_yld = [sin(Az_yld)*cos(El_yld); cos(Az_yld)*cos(El_yld); sin(El_yld)];
dir_chj = [sin(Az_chj)*cos(El_chj); cos(Az_chj)*cos(El_chj); sin(El_chj)];

