function [x, y, z] = cramer_rule(A, p)
    % A 是系数矩阵，p 是结果向量
    % 假设 A 为 3x3 矩阵，p 为 3x1 向量

    % 计算矩阵 A 的行列式
    det_A = det(A);
    
    if det_A == 0
        error('矩阵 A 的行列式为零，无法使用克莱姆法则求解');
    end
    
    % 构造矩阵 A_x, A_y, A_z
    A_x = A; A_y = A; A_z = A;
    A_x(:,1) = p;  % 替换 A 的第一列
    A_y(:,2) = p;  % 替换 A 的第二列
    A_z(:,3) = p;  % 替换 A 的第三列
    
    % 计算 A_x, A_y, A_z 的行列式
    det_A_x = det(A_x);
    det_A_y = det(A_y);
    det_A_z = det(A_z);
    
    % 使用克莱姆法则计算 x, y, z
    x = det_A_x / det_A;
    y = det_A_y / det_A;
    z = det_A_z / det_A;
end
