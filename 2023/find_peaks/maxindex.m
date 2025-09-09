function max_index = maxindex(vector)
    % 提取实部部分
    
    max_value = max(vector);
    % 找到最大值对应的索引
    max_index = find(vector == max_value);
end