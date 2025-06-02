function T = KalmanFilter(s,Q,R)
% --------------------------------
% T = KalmanFilter(s,Q,R)
% s 位移，列向量
% Q 状态误差协方差矩
% R 观测误差协方差矩
% --------------------------------
% 建立质点运动方程:
% S(k+1) = 1*S(k) + T*V(k) + 0.5*T^2*a
% V(k+1) = 0*S(k) + 1*V(k) + T*a
% 建立观测方程:
% y(k+1) = S(k+1) + v(k+1)
% 即：
% X(k+1) = A * X(k) + G*w(k+1); 预测模型
% y(k+1) = H * X(k+1) + v(k+1); 观测模型
N = length(s);
T = 1; % 采样间隔默认为1
A = [1 T;0 1]; % 状态转移矩阵
G = [T^2/2;T]; % 控制量矩阵
H = [1 0]; % 观测矩阵
% 初始化第一个状态
Xu = [s(1); 0];
Pu = [0 0;0 0];
I = [1 0;0 1];
T = zeros(N,1);
for i = 2:N
Xp = A * Xu;
Pp = A * Pu * A' + G * Q * G';
K = Pp * H' * ( H * Pp * H' + R)^-1;
Xu = ( I - K * H ) * Xp + K * s(i);
Pu = ( I - K * H ) * Pp;
T(i) = Xu(1);
end
end