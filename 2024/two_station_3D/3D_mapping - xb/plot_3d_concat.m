% 加载数据
old_method_match_result = load('match_results_all_4096.mat'); 
old_method_S_result = load("S_results_all_4096.mat");
new_method_match_result = load("all_match_result_0.8_0.1_0.1.mat"); 
new_method_S_result = load("all_S_result_0.8_0.1_0.1.mat");

% 定义范围
x_range = [-10000, 8000];
y_range = [-10000, 0];
z_range = [0, 10000];

%% --- 筛选 new_method 数据 ---
% 筛选条件
conditions_new_method = ([new_method_match_result.all_match_results.chi_square_red] < 1000) & ...
            ([new_method_match_result.all_match_results.dlta] < 40000) & ...
             ([new_method_match_result.all_match_results.yld_start_loc] > 3e8) & ...
             ([new_method_match_result.all_match_results.yld_start_loc] < 4.5e8) & ...
             ([new_method_match_result.all_match_results.r_gccs] > 0.2) & ...
             (abs([new_method_match_result.all_match_results.R3_value]) < 800);

% 应用筛选
filtered_match_indices_new_method = find(conditions_new_method);
filtered_S_temp_new_method = new_method_S_result.all_S_results(filtered_match_indices_new_method, :);
filtered_match_result_temp_new_method = new_method_match_result.all_match_results(filtered_match_indices_new_method); 

% 进一步筛选空间范围
range_condition_s_new_method = filtered_S_temp_new_method(:,1) >= x_range(1) & filtered_S_temp_new_method(:,1) <= x_range(2) & ...
    filtered_S_temp_new_method(:,2) >= y_range(1) & filtered_S_temp_new_method(:,2) <= y_range(2) & ...
    filtered_S_temp_new_method(:,3) >= z_range(1) & filtered_S_temp_new_method(:,3) <= z_range(2);

filtered_S_new_method = filtered_S_temp_new_method(range_condition_s_new_method, :);
filtered_match_result_new_method = filtered_match_result_temp_new_method(range_condition_s_new_method);

% 检查 yld_start_loc 是否可用
if ~isempty(filtered_match_result_new_method) && isfield(filtered_match_result_new_method, 'yld_start_loc') && isnumeric([filtered_match_result_new_method.yld_start_loc])
    time_colors_new_method = [filtered_match_result_new_method.yld_start_loc]';
else
    disp('警告: filtered_match_result_new_method 为空，或 yld_start_loc 不可用/非数值类型。将按索引着色。');
    time_colors_new_method = (1:size(filtered_S_new_method,1))';
end

%% --- 筛选 old_method 数据 ---
% 筛选条件
conditions_old_method = ([old_method_match_result.all_match_results.dlta] < 50000) & ...
             ([old_method_match_result.all_match_results.yld_start_loc] > 4.5e8) & ...
             ([old_method_match_result.all_match_results.r_gccs] > 0.1);

% 应用筛选
filtered_match_indices_old_method = find(conditions_old_method);
filtered_S_temp_old_method = old_method_S_result.all_S_results(filtered_match_indices_old_method, :);
filtered_match_result_temp_old_method = old_method_match_result.all_match_results(filtered_match_indices_old_method); 

% 进一步筛选空间范围
range_condition_s_old_method = filtered_S_temp_old_method(:,1) >= x_range(1) & filtered_S_temp_old_method(:,1) <= x_range(2) & ...
    filtered_S_temp_old_method(:,2) >= y_range(1) & filtered_S_temp_old_method(:,2) <= y_range(2) & ...
    filtered_S_temp_old_method(:,3) >= z_range(1) & filtered_S_temp_old_method(:,3) <= z_range(2);

filtered_S_old_method = filtered_S_temp_old_method(range_condition_s_old_method, :);
filtered_match_result_old_method = filtered_match_result_temp_old_method(range_condition_s_old_method);

% 检查 yld_start_loc 是否可用
if ~isempty(filtered_match_result_old_method) && isfield(filtered_match_result_old_method, 'yld_start_loc') && isnumeric([filtered_match_result_old_method.yld_start_loc])
    time_colors_old_method = [filtered_match_result_old_method.yld_start_loc]';
else
    disp('警告: filtered_match_result_old_method 为空，或 yld_start_loc 不可用/非数值类型。将按索引着色。');
    time_colors_old_method = (1:size(filtered_S_old_method,1))';
end

%% --- 合并 new_method 和 old_method 数据 ---
filtered_S = [filtered_S_new_method; filtered_S_old_method];  % new_method 在前，old_method 在后
time_colors = [time_colors_new_method; time_colors_old_method];  % 合并颜色数据

marker_size = 1;


x = filtered_S(:, 1);
y = filtered_S(:, 2); 
z = filtered_S(:, 3); 

% --- 4.1 三维散点图 (按时间着色) ---
figure;
%绘制三维散点图：x, y, z 为坐标，
scatter3(x, y, z, marker_size, time_colors, 'filled');
xlabel('X (东)'); 
ylabel('Y (北)'); 
zlabel('Z (上)'); 
title('源的三维空间分布');
xlim(x_range);
ylim(y_range);
zlim(z_range);
grid on; 
axis equal;
colorbar; 
colormap(gca, 'hsv');
daspect([1 1 1]); 

%% --- 5. 自定义数据显示功能 ---
% 获取当前图窗的数据光标模式对象
dcm = datacursormode(gcf); 
% 启用数据光标模式
datacursormode on;

% 设置自定义的更新函数
% @(obj, event_obj) 是一个匿名函数，它会在数据光标事件发生时执行
% event_obj 包含了事件信息，其中 event_obj.DataIndex 是被点击数据点的索引
dcm.UpdateFcn = @(obj, event_obj) ...
    {sprintf('X (东): %.2f', event_obj.Position(1)), ...
     sprintf('Y (北): %.2f', event_obj.Position(2)), ...
     sprintf('Z (上): %.2f', event_obj.Position(3)), ...
     sprintf('yld-start-loc: %.0f', time_colors(event_obj.DataIndex))};

