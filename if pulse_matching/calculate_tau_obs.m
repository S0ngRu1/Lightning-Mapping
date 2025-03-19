% 定义计算τij的理想值τ_ij^obs的函数
function tau_ij_obs = calculate_tau_obs(cos_alpha, cos_beta)
    %从化局
    angle12 = -2.8381;
    angle13 = 28.2006;
    angle23 = 87.3358;
    d12 = 41.6496;
    d13 = 48.5209;
    d23 = 25.0182;
    % 使用式(3)计算τij的理想值τ_ij^obs
    tau_ij_obs(1) = (cos_alpha * sind(angle12) + cos_beta * cosd(angle12)) * d12 / 0.299792458;
    tau_ij_obs(2) = (cos_alpha * sind(angle13) + cos_beta * cosd(angle13)) * d13 / 0.299792458;
    tau_ij_obs(3) = (cos_alpha * sind(angle23) + cos_beta * cosd(angle23)) * d23 / 0.299792458;
%     %引雷场
%     angle12 = -110.8477;
%     angle13 = -65.2405;
%     angle23 = -19.6541;
%     d12 = 24.9586;
%     d13 = 34.9335;
%     d23 = 24.9675;
%     % 使用式(3)计算τij的理想值τ_ij^obs
%     tau_ij_obs(1) = (cos_alpha * sind(angle12) + cos_beta * cosd(angle12)) * d12 / 0.299792458;
%     tau_ij_obs(2) = (cos_alpha * sind(angle13) + cos_beta * cosd(angle13)) * d13 / 0.299792458;
%     tau_ij_obs(3) = (cos_alpha * sind(angle23) + cos_beta * cosd(angle23)) * d23 / 0.299792458;

end