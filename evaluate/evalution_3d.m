%% 1. 加载数据
% 假设你已经手动筛选出了某一段连续通道的数据，保存为 subset_data.csv
% 或者从总表中通过索引提取
data = readtable('..\2024\two_station_3D\results\3d_win512_cost_cal_yld_chj_dtoa_3.6e8_5.6e8_with_initial_results.csv');

target_start_loc = 4.68e8;  % 起始采样点位置 (需要根据实际数据修改)
target_end_loc   = 4.712e8;  % 结束采样点位置 (需要根据实际数据修改)

% 使用逻辑索引进行筛选
segment_idx = data.yld_start_loc >= target_start_loc & data.yld_start_loc <= target_end_loc;
segment = data(segment_idx, :);

% 提取坐标
t = segment.yld_start_loc; 
X_tri = segment.x_tri; Y_tri = segment.y_tri; Z_tri = segment.z_tri;
X_dtoa = segment.x_dtoa; Y_dtoa = segment.y_dtoa; Z_dtoa = segment.z_dtoa;


%% 3. 计算并对比平滑度 (Scattering)
poly_order = 4; % 论文设定为 5 阶

% --- 三角测量 (Triangulation) ---
[std_x_tri, res_x_tri] = evaluate_smoothness(t, X_tri, poly_order);
[std_y_tri, res_y_tri] = evaluate_smoothness(t, Y_tri, poly_order);
[std_z_tri, res_z_tri] = evaluate_smoothness(t, Z_tri, poly_order);

% --- DTOA 优化 (Refined) ---
[std_x_dtoa, res_x_dtoa] = evaluate_smoothness(t, X_dtoa, poly_order);
[std_y_dtoa, res_y_dtoa] = evaluate_smoothness(t, Y_dtoa, poly_order);
[std_z_dtoa, res_z_dtoa] = evaluate_smoothness(t, Z_dtoa, poly_order);

%% 4. 输出评估报告
fprintf('=== 通道平滑度评估结果 (1-sigma Scattering) ===\n');
fprintf('单位: 米 (m)\n\n');

fprintf('方向\t 三角测量(Tri)\t DTOA优化\t 提升倍数\n');
fprintf('Easting(X)\t %.2f\t\t %.2f\t\t %.1fx\n', std_x_tri, std_x_dtoa, std_x_tri/std_x_dtoa);
fprintf('Northing(Y)\t %.2f\t\t %.2f\t\t %.1fx\n', std_y_tri, std_y_dtoa, std_y_tri/std_y_dtoa);
fprintf('Altitude(Z)\t %.2f\t\t %.2f\t\t %.1fx\n', std_z_tri, std_z_dtoa, std_z_tri/std_z_dtoa);

%% 5. 绘图对比 (复现论文图7)
figure('Position', [100, 100, 1000, 400]);
subplot(1,2,1);
scatter(X_tri, Y_tri, 10, 'filled', 'MarkerFaceColor', 'b'); hold on;
plot(polyval(polyfit(t, X_tri, 5), t), polyval(polyfit(t, Y_tri, 5), t), 'k', 'LineWidth', 1.5);
title(['Triangulation (Scatter: ' sprintf('%.1fm', mean([std_x_tri, std_y_tri])) ')']);
xlabel('Easting (m)'); ylabel('Northing (m)'); grid on; axis equal;

subplot(1,2,2);
scatter(X_dtoa, Y_dtoa, 10, 'filled', 'MarkerFaceColor', 'r'); hold on;
plot(polyval(polyfit(t, X_dtoa, 5), t), polyval(polyfit(t, Y_dtoa, 5), t), 'k', 'LineWidth', 1.5);
title(['DTOA Optimized (Scatter: ' sprintf('%.1fm', mean([std_x_dtoa, std_y_dtoa])) ')']);
xlabel('Easting (m)'); ylabel('Northing (m)'); grid on; axis equal;



%% 2. 定义评估函数
function [std_res, residuals] = evaluate_smoothness(t, pos, poly_order)
    % pos: 坐标数据向量 (如 X)
    % poly_order: 拟合阶数，论文中使用 5 阶
    
    % 剔除 NaN
    valid_idx = ~isnan(pos);
    t_valid = t(valid_idx);
    pos_valid = pos(valid_idx);
    
    if length(pos_valid) < poly_order + 2
        std_res = NaN; residuals = []; return;
    end

    % 多项式拟合
    p = polyfit(t_valid, pos_valid, poly_order);
    pos_fit = polyval(p, t_valid);
    
    % 计算残差
    residuals = pos_valid - pos_fit;
    
    % 计算标准差 (RMSE)
    std_res = std(residuals);
end
