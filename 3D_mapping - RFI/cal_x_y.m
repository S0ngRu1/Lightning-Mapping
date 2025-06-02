% 输入两地的经纬度（单位：度）
lat1= 23.6392983;%引雷点的纬度
lon1= 113.5957504;%引雷点的经度
lat2 = 23.5684654;
lon2 = 113.6153785;
 
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



% % 输入经纬度坐标（WGS84坐标系）
% CHJ_lon = 113.6153785;
% CHJ_lat = 23.5684654;
% YLD_lon = 113.5957504;
% YLD_lat = 23.6392983;
% 
% CHJ_lat = 23.5684654;
% CHJ_lon = 113.6153785;
% YLD_lat = 23.6392983;
% YLD_lon = 113.5957504;
% 
% % 地球半径近似值 (单位：米)
% R = 6371000; % 地球半径，单位：米
% 
% % 计算两点的经纬度差值
% delta_lat = CHJ_lat - YLD_lat;
% delta_lon = CHJ_lon - YLD_lon;
% 
% % 将经纬度差值转换为米（单位：米）
% lat_distance = delta_lat * 111000; % 纬度差距转换为米
% lon_distance = delta_lon * 111000 * cosd(YLD_lat); % 经度差距转换为米
% 
% % 以YLD为原点，将CHJ坐标变换为新坐标系
% X_CHJ = lon_distance; % 东西方向 (X轴)
% Y_CHJ = lat_distance; % 南北方向 (Y轴)
% 
% % 输出结果
% fprintf('CHJ在以YLD为原点的坐标系中的位置：X = %.2f 米, Y = %.2f 米\n', X_CHJ, Y_CHJ);



