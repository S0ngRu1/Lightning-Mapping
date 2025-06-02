function valid_indices = find_matching_indices(valid_chj_start_loc, sample_point_gap, yld_start_loc)
    % 计算有效区间的上下限
    lower_bound = valid_chj_start_loc - sample_point_gap;
    upper_bound = valid_chj_start_loc + sample_point_gap;

    % 查找区间 [lower_bound, upper_bound] 在 yld_start_loc 中的位置范围
    valid_indices = find_in_range(yld_start_loc, lower_bound, upper_bound);
end

function valid_indices = find_in_range(yld_start_loc, lower_bound, upper_bound)
    % 利用二分查找找到区间的左右边界
    
    % 查找第一个大于等于 lower_bound 的元素索引
    left_index = binary_search(yld_start_loc, lower_bound, 'left');
    
    % 查找第一个大于 upper_bound 的元素索引
    right_index = binary_search(yld_start_loc, upper_bound, 'right');
    
    % 获取区间内所有的有效索引
    valid_indices = left_index:right_index-1;
end

function index = binary_search(arr, target, direction)
    % 二分查找：找大于或小于目标值的索引
    
    left = 1;
    right = length(arr);
    index = -1;
    
    while left <= right
        mid = floor((left + right) / 2);
        
        if arr(mid) < target
            left = mid + 1;
        elseif arr(mid) > target
            right = mid - 1;
        else
            index = mid;
            if strcmp(direction, 'left')
                right = mid - 1;
            elseif strcmp(direction, 'right')
                left = mid + 1;
            end
        end
    end
    
    if strcmp(direction, 'left')
        index = left;
    elseif strcmp(direction, 'right')
        index = right + 1;
    end
end
