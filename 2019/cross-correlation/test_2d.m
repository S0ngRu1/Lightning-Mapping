%读取数据
signal_length = 1e6;
ch1 = read_signal('20190604164852.7960CH1.dat',signal_length);
ch2 = read_signal('20190604164852.7960CH2.dat',signal_length);

window_length = 1024;
windows =1:256:signal_length-window_length+1;
d = 20;
c = 299.552816;

fileID = fopen('result1_2d.txt', 'w');
% 写入第一行的数据介绍
fprintf(fileID, '%-13s%-15s%-15s%-15s%-15s\n', ...
    'Start_loc','tau', 'cos_alpha', 'Elevation', 'corr12');

%遍历所有窗口
for  wi = 1:numel(windows)
%   读信号
    signal1 = ch1(windows(wi):windows(wi)+window_length-1);
    signal2 = ch2(windows(wi):windows(wi)+window_length-1);

% 
    filtered_signal1 = filter_fft(signal1,40 ,80);
    filtered_signal2 = filter_fft(signal2,40 ,80);


    % 对滤波后的信号应用窗函数
    windowed_signal1 = windowsignal(filtered_signal1);
    windowed_signal2 = windowsignal(filtered_signal2);

    %处理后的信号
    ch1_new = windowed_signal1;
    ch2_new = windowed_signal2;
    [tau,R,lag] = gccphat(ch1_new,ch2_new,200);
    r = max(R);
    cos_alpha = tau*c/d;
    if cos_alpha > 1 || cos_alpha <0
        continue;
    end
    El = acos(cos_alpha);
    % 将弧度转换为角度
    El_deg = rad2deg(El);
    % 写入计算后的数据
        fprintf(fileID, '%-13d%-15.6f%-15.6f%-15.6f%-15.6f\n', ...
             400000000+windows(wi), tau, cos_alpha,El_deg, r);

end
% 关闭文件
fclose(fileID);




