% % 定义目标函数
% function F = objective(x,t12,t13,t23,type)
%     % 提取待优化的变量
%     cos_alpha = x(1);
%     cos_beta = x(2);
% 
%     % 计算τij的理想值τ_ij^obs
%     tau_ij_obs = calculate_tau_obs(cos_alpha, cos_beta,type);
%     % 计算Δt12, Δt13, Δt23
%     delta_t12 = delta_t(t12,tau_ij_obs(1));
%     delta_t13 = delta_t(t13,tau_ij_obs(2));
%     delta_t23 = delta_t(t23,tau_ij_obs(3));
% 
%     % 计算目标函数，即式(4)
%     F = (delta_t12^2 + delta_t13^2 + delta_t23^2) / 75;
% end

% 定义目标函数 (正确版本)
function F = objective(x, t12_meas, t13_meas, t23_meas, type)
    % 提取待优化的变量
    cos_alpha = x(1);
    cos_beta = x(2);

    % 计算τij的理论值 τ_model (我将 obs 改为 model，语义更清晰)
    tau_model = calculate_tau_obs(cos_alpha, cos_beta, type);
    
    % t12, t13, t23 是测量的时延 (measurement)
    % tau_model(1), tau_model(2), tau_model(3) 是根据当前 x 计算出的理论时延
    
    % 计算残差向量
    residual12 = t12_meas - tau_model(1);
    residual13 = t13_meas - tau_model(2);
    residual23 = t23_meas - tau_model(3);

    % 返回残差向量 F
    % lsqnonlin 会自动最小化 sum(F.^2)
    F = [residual12; residual13; residual23];
end