%函数：对主窗口进行上采样
function new_signal = upsampling(original_signal,upsampling_factor)

    % 原信号
    original_x = (1:numel(original_signal))';
    original_y = original_signal;
    % 上采样后的采样点数
    upsampled_length = length(original_x) * upsampling_factor;
    % 上采样后的采样点的 x 坐标
    upsampled_x = linspace(1, length(original_x), upsampled_length);
    % 使用多项式插值对原信号进行上采样
    interpolated_signal = interp1(original_x, original_y, upsampled_x, 'spline');
    new_signal = [upsampled_x; interpolated_signal];
end

% function new_signal = upsampling(original_signal, upsampling_factor, method)
    % 输入参数：
    % original_signal - 原始信号
    % upsampling_factor - 上采样因子
    % method - 插值方法，可选值：'linear', 'spline', 'polyfit'
    
%     % 原信号的采样点位置
%     original_x = (1:numel(original_signal))';
%     original_y = original_signal;
% 
%     % 计算上采样后的信号长度和位置
%     upsampled_length = length(original_x) * upsampling_factor;
%     upsampled_x = linspace(1, length(original_x), upsampled_length)';
% 
%     switch method
%         case 'linear'
%             % 线性插值
%             interpolated_signal = interp1(original_x, original_y, upsampled_x, 'linear');
%             
%         case 'spline'
%             % 样条插值
%             interpolated_signal = interp1(original_x, original_y, upsampled_x, 'spline');
%             
%         case 'polyfit'
%             % 多项式拟合插值
%             poly_order = 1; % 设定多项式的阶数
%             p = polyfit(original_x, original_y, poly_order);
%             interpolated_signal = polyval(p, upsampled_x);
%             
%         otherwise
%             error('Unsupported interpolation method. Use ''linear'', ''spline'', or ''polyfit''.');
%     end
% 
%     % 返回上采样后的信号
%     new_signal = [upsampled_x, interpolated_signal];
% end
% 
