signal_length = 4e6; 
start_read_loc_yld = 3.92e8;
start_read_loc_chj = start_read_loc_yld+34151156;
yld_signal = read_signal('../20240822165932.6610CH1.dat',signal_length,start_read_loc_yld);
chj_signal = read_signal('../2024 822 85933.651462CH1.dat',signal_length,start_read_loc_chj);
x = downsample(yld_signal,1);
y = downsample(chj_signal,1);
% plot(x);
subplot(2,1,1);plot(x);title('yld');xlabel('采样点');ylabel('振幅');
subplot(2,1,2);plot(y);title('chj');xlabel('采样点');ylabel('振幅');
figure;
plot(chj_signal);title('chj');xlabel('采样点');ylabel('振幅');
filtered_chj_signal1 = filter_bp(chj_signal,30e6,80e6,5);
plot(filtered_chj_signal1);title('chj');xlabel('采样点');ylabel('振幅');


signal_length = 12820000; 
index = 790;
yld_signal = read_signal('../20240822165932.6610CH1.dat',signal_length,all_start_signal_loc_yld(index));
chj_signal = read_signal('../2024 822 85933.651462CH1.dat',signal_length,all_start_signal_loc_chj(index));
filtered_yld_signal1 = filter_bp(yld_signal,30e6,80e6,5);
filtered_chj_signal2 = filter_bp(chj_signal,30e6,80e6,5);
subplot(2,1,1);plot(filtered_yld_signal1);title('yld');xlabel('采样点');ylabel('振幅');
subplot(2,1,2);plot(filtered_chj_signal2);title('chj');xlabel('采样点');ylabel('振幅');


% 定义系数矩阵 A
A = [4246610, 1;
     (25640000+10163123), 1];

% 定义常数项列向量 B
B = [1790;
     9300];

% 使用左除求解
solutions = A \ B;

% 提取解
x_sol = solutions(1);
b_sol = solutions(2);

% 显示结果
fprintf('使用矩阵除法求解:\n');
fprintf('x = %f\n', x_sol);
fprintf('b = %f\n', b_sol);
fprintf('x (科学计数法) = %e\n', x_sol);
fprintf('b (科学计数法) = %e\n', b_sol);

% 验证解 (可选)
error1 = A(1,1)*x_sol + A(1,2)*b_sol - B(1);
error2 = A(2,1)*x_sol + A(2,2)*b_sol - B(2);
fprintf('方程1的误差: %e\n', error1);
fprintf('方程2的误差: %e\n', error2);



% 求两个点的距离
x1 = 0;
y1 = 0;
x2 = -2.0622;
y2 = 41.5985;
d = sqrt((x1-x2)^2+(y1-y2)^2);
