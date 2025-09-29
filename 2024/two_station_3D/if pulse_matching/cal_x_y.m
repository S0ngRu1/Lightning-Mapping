% 输入两地的经纬度（单位：度）
lat1 = 23.6392983; % 第一个点的纬度
lon1 = 113.5957504; % 第一个点的经度
lat2 = 23.5684654; % 第二个点的纬度
lon2 = 113.6153785; % 第二个点的经度

% 将经纬度转换为弧度
lat1_rad = deg2rad(lat1);
lon1_rad = deg2rad(lon1);
lat2_rad = deg2rad(lat2);
lon2_rad = deg2rad(lon2);

% 地球半径（单位：米）
R = 6371000;

% 计算第二点相对于第一个点的平面直角坐标
d_lat = lat2_rad - lat1_rad;
d_lon = lon2_rad - lon1_rad;
x = R * d_lon * cos((lat1_rad + lat2_rad) / 2); % 东西方向
y = R * d_lat; % 南北方向

% 计算相对于 y 轴的夹角（顺时针方向）
theta_rad = atan2(x, y); % 弧度制，相对 y 轴角度
theta_deg = rad2deg(theta_rad); % 转换为角度制

% 将角度调整为顺时针从 y 轴开始
if theta_deg < 0
    theta_deg = theta_deg + 360;
end

% 输出结果
disp(['第二点相对于第一个点的坐标为：x = ', num2str(x), ' 米, y = ', num2str(y), ' 米']);
disp(['相对于 y 轴的夹角为（顺时针）：', num2str(theta_deg), ' 度']);

