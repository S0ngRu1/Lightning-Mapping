function max_index = maxindex(vector)
    % 提取实部部分
    real_vector = real(vector);
    % 找到实部大于零的元素
    positive_values = real_vector(real_vector > 0);
    % 找到实部大于零的元素中的最大值
    max_value = max(positive_values);
    % 找到最大值对应的索引
    max_index = find(real_vector == max_value);
end