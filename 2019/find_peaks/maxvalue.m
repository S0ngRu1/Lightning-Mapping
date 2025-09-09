function max_value = maxvalue(vector)
    % 提取实部部分
    real_vector = real(vector);
    % 找到实部大于零的元素
    positive_values = real_vector(real_vector > 0);
    % 找到实部大于零的元素中的最大值
    max_value = max(positive_values);
end