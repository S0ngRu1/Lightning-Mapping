%对互相关函数进行上采样
function upsampling_gcc = upsampling_gc(r,lag,upsampling_factor)
    % 创建插值函数对象               
    interpolant = griddedInterpolant(lag, r, 'spline');
    % 定义新的插值网格
    new_x = linspace(-numel(r)/2, numel(r)/2, numel(r)*upsampling_factor);
    % 执行插值上采样
    interpolated_signal = interpolant(new_x);
    upsampling_gcc = [new_x; interpolated_signal]';
end