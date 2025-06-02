function [start_loc, azimuth, elevation, Rcorr, t123]  = get_2d_result_single_window(start_loc,ch1,ch2,ch3 )
% 从化局
N = 3;
d12 = 41.6496;
d13 = 48.5209;
d23 = 25.0182;
angle12 = -2.8381;
angle13 = 28.2006;
angle23 = 87.3358;
c = 0.299792458;
fs = 200e6;

[ch1_up, ch2_up, ch3_up] = deal(...
    upsampling(ch1, 50)', ...
    upsampling(ch2, 50)', ...
    upsampling(ch3, 50)');
ch1_upsp = ch1_up(:,2);
ch2_upsp = ch2_up(:,2);
ch3_upsp = ch3_up(:,2);


[r12_gcc,lags12_gcc] = xcorr(ch1_upsp,ch2_upsp,'normalized');
[r13_gcc,lags13_gcc] = xcorr(ch1_upsp,ch3_upsp,'normalized');
[r23_gcc,lags23_gcc] = xcorr(ch2_upsp,ch3_upsp,'normalized');
R12_gcc = max(r12_gcc);
R13_gcc = max(r13_gcc);
R23_gcc = max(r23_gcc);
t12_gcc = cal_tau(r12_gcc,lags12_gcc');
t13_gcc = cal_tau(r13_gcc,lags13_gcc');
t23_gcc = cal_tau(r23_gcc,lags23_gcc');


t12 = t12_gcc *0.1;
t13 = t13_gcc *0.1+2;
t23 = t23_gcc *0.1+2;
cos_beta_0 =((c*t13*d12*sind(angle12))-(c*t12*sind(angle13)*d13))/(d13*d12*sind(angle12-angle13)) ;
cos_alpha_0 = ((c*t12)/d12-cos_beta_0*cosd(angle12))/sind(angle12);
if abs(cos_beta_0)>1 || abs(cos_alpha_0)>1
    start_loc = 0;
    t123 = 0;
    azimuth = 0;
    elevation = 0;
    Rcorr = 0;
    return
end
x0 = [cos_alpha_0,cos_beta_0];
% 调用lsqnonlin函数进行优化
options = optimoptions('lsqnonlin', 'MaxIter', 1000, 'TolFun', 1e-6);
x = lsqnonlin(@(x) objective(x, t12, t13, t23), x0, [-1 -1], [1 1], options);
% 输出最优的cos(α)和cos(β)值
cos_alpha_opt = x(1);
cos_beta_opt = x(2);
if abs(cos_alpha_opt)>1 || abs(cos_beta_opt)>1
    start_loc = 0;
    t123 = 0;
    azimuth = 0;
    elevation = 0;
    Rcorr = 0;
    return
end
Az = atan2( cos_alpha_opt,cos_beta_opt);
if abs(cos_beta_opt/cos(Az)) > 1
    start_loc = 0;
    t123 = 0;
    azimuth = 0;
    elevation = 0;
    Rcorr = 0;

    return
end
El = acos( cos_beta_opt/cos(Az) );
% 将弧度转换为角度
azimuth = rad2deg(Az);
elevation = rad2deg(El);
if azimuth < 0
    azimuth = azimuth + 360;
end

t123 = t12 + t23 - t13;
Rcorr = (R12_gcc + R13_gcc + R23_gcc) / 3;
if t123>1 || Rcorr<0.3
    start_loc = 0;
    t123 = 0;
    azimuth = 0;
    elevation = 0;
    Rcorr = 0;
end

end