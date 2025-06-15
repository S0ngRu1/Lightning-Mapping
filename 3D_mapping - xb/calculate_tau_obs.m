% 定义计算τij的理想值τ_ij^obs的函数
function tau_ij_obs = calculate_tau_obs(cos_alpha, cos_beta, type)
    % 初始化输出变量
    tau_ij_obs = zeros(1, 3);

    % 根据 type 参数选择不同的参数集
    if strcmp(type, 'chj') % 从化局
        angle12 = -2.8381;
        angle13 = 50.3964;
        angle23 = 120.6568;
        d12 = 41.6496;
        d13 = 36.9015;
        d23 = 35.4481;
    elseif strcmp(type, 'yld') % 引雷场
        angle12 = -110.8477;
        angle13 = -65.2405;
        angle23 = -19.6541;
        d12 = 24.9586;
        d13 = 34.9335;
        d23 = 24.9675;
    else
        error('未知的类型：%s', type);
    end

    % 使用式(3)计算τij的理想值τ_ij^obs
    tau_ij_obs(1) = (cos_alpha * sind(angle12) + cos_beta * cosd(angle12)) * d12 / 0.299792458;
    tau_ij_obs(2) = (cos_alpha * sind(angle13) + cos_beta * cosd(angle13)) * d13 / 0.299792458;
    tau_ij_obs(3) = (cos_alpha * sind(angle23) + cos_beta * cosd(angle23)) * d23 / 0.299792458;
end