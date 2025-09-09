%% 使用EEMD对信号进行降噪和滤波
function filtered_signal = filter_eemd(sig)
    %EEMD分解
    Nstd = 0.2; %Nstd为附加噪声标准差与Y标准差之比
    NE = 20;   %NE为对信号的平均次数
    imf = eemd(sig,Nstd,NE);
    %信号重构
    indices = [2,5,6,7,8,9,10,11];  
    filtered_imfs = imf(:,indices);
    filtered_signal = sum(filtered_imfs, 2);
end