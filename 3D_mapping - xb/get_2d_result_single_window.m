function [start_loc,t12,t13,t23, Az_deg, El_deg, Rcorr, t123]  = get_2d_result_single_window(start_loc,ch1,ch2,ch3 ,type)
if strcmp(type, 'chj')
    % 从化局
    d12 = 41.6496;
    d13 = 36.9015;
    d23 = 35.4481;
    angle12 = -2.8381;
    angle13 = 50.3964;
    angle23 = 120.6568;
elseif strcmp(type, 'yld')
    % 引雷点
    angle12 = -110.8477;
    angle13 = -65.2405;
    angle23 = -19.6541;
    d12 = 24.9586;
    d13 = 34.9335;
    d23 = 24.9675;

end
start_loc = 0;
Az_deg = 0;
El_deg = 0;
Rcorr = 0;
t123 = 100;
c = 0.299792458;
[ch1_up, ch2_up, ch3_up] = deal(...
    upsampling(ch1, 50)', ...
    upsampling(ch2, 50)', ...
    upsampling(ch3, 50)');
ch1_upsp = ch1_up(:,2);
ch2_upsp = ch2_up(:,2);
ch3_upsp = ch3_up(:,2);
%互相关
[r12_gcc,lags12_gcc] = xcorr(ch1_upsp,ch2_upsp,'normalized');
[r13_gcc,lags13_gcc] = xcorr(ch1_upsp,ch3_upsp,'normalized');
[r23_gcc,lags23_gcc] = xcorr(ch2_upsp,ch3_upsp,'normalized');
R12_gcc = max(r12_gcc);
R13_gcc = max(r13_gcc);
R23_gcc = max(r23_gcc);
t12_gcc = cal_tau(r12_gcc,lags12_gcc');
t13_gcc = cal_tau(r13_gcc,lags13_gcc');
t23_gcc = cal_tau(r23_gcc,lags23_gcc');
%                 从化局
t12 = t12_gcc *0.1;
t13 = t13_gcc *0.1+1.600061;
t23 = t23_gcc *0.1+1.600061;
cos_beta_0 =((c*t13*d12*sind(angle12))-(c*t12*sind(angle13)*d13))/(d13*d12*sind(angle12-angle13)) ;
cos_alpha_0 = ((c*t12)/d12-cos_beta_0*cosd(angle12))/sind(angle12);
if abs(cos_beta_0)>1 || abs(cos_alpha_0)>1
    return;
end
x0 = [cos_alpha_0,cos_beta_0];
% 调用lsqnonlin函数进行优化
options = optimoptions('lsqnonlin', 'MaxIter', 1000, 'TolFun', 1e-6);
x = lsqnonlin(@(x) objective(x, t12, t13, t23,type), x0, [-1 -1],[1 1], options);
% 输出最优的cos(α)和cos(β)值
cos_alpha_opt = x(1);
cos_beta_opt = x(2);
if abs(cos_alpha_opt)>1 || abs(cos_beta_opt)>1
    return;
end
Az = atan2( cos_alpha_opt,cos_beta_opt);
if abs(cos_beta_opt/cos(Az)) > 1
    return;
end
El = acos( cos_beta_opt/cos(Az) );
% 将弧度转换为角度
Az_deg = rad2deg(Az);
El_deg = rad2deg(El);
if Az_deg < 0
    Az_deg = Az_deg + 360;
end

t123 = t12 + t23 - t13;
Rcorr = (R12_gcc + R13_gcc + R23_gcc) / 3;
end