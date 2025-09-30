function shifted_signal = shift_signal(signal, shift_amount)

    % 使用 circshift 进行平移
    shifted_signal = circshift(signal, shift_amount);
    % 如果是向左平移，右侧补零；如果是向右平移，左侧补零
    if shift_amount < 0
        shifted_signal(end+shift_amount+1:end) = 0;
    else
        shifted_signal(1:shift_amount) = 0;
    end
    
end
