% function tau = cal_tau(R, lag)
%     % 从数据中找到y的最大值及其索引
%     [~, max_index] = max(R);
%     tau = lag(max_index,1);
% end
% 
function tau = cal_tau(R, lag) % 新的、高精度的版本
    [~, max_idx] = max(R);

    % 确保峰值不在数组的边缘，否则无法取到3个点
    if max_idx == 1 || max_idx == length(R)
        tau = lag(max_idx);
        return;
    end

    % 提取峰值点 (y2) 和它左右相邻的两个点 (y1, y3)
    y1 = R(max_idx - 1);
    y2 = R(max_idx);
    y3 = R(max_idx + 1);

    % 抛物线顶点横坐标的偏移量公式： p = (y1 - y3) / (2 * (y1 - 2*y2 + y3))
    % p 是相对于中心点 max_idx 的亚采样偏移量
    % 注意：要处理分母为0或非常小的情况，避免计算错误
    denominator = 2 * (y1 - 2*y2 + y3);
    if abs(denominator) < 1e-9
        p = 0; % 如果分母太小（例如，平顶），则不进行偏移
    else
        p = (y1 - y3) / denominator;
    end
    
    % 计算最终的精确时延
    % lag是等差数列，可以直接用 p 乘以步长
    time_step = lag(2) - lag(1);
    tau = lag(max_idx) + p * time_step;
end