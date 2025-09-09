%%
% 读取原始信号
signal_length = 2e3;
sig = read_signal('20190604164852.7960CH3.dat',signal_length);
fs = 200;
t = (0:signal_length-1) * 1/fs; % 将时间单位从秒改为微秒
figure('color','white')
subplot(211)
plot(t, sig, 'k') % 绘制原始信号
xlabel('Time (\mu s)') % 设置 x 轴标签单位为微秒
ylabel('Amplitude')
title('Original Signal')

%%
% 傅里叶变换，画出频谱图

y=fft(sig);%傅里叶变换得到一个复数
n = length(y);
x = (0:n/2-1) * (fs/n);
Ay=abs(y);%取模
Ayy=Ay*2/signal_length;%转换成实际的幅值
figure(2)
plot(x,Ayy(1:signal_length/2))
xlabel('Frequency (MHz)')
ylabel('Amplitude')

%%
%设置两个频率参数 f1 和 f2 分别为 8 和 15；
f1=40;
f2=80;
%创建一个长度与输入信号 y 相同的零向量 yy；
yy=zeros(1,length(y));
% 使用 for 循环遍历信号 y 的每个采样点（m 表示当前的采样点索引，从0到 N-1）；
for m=1:n-1
%     判断当前采样点对应的频率是否在 8Hz 到 15Hz 范围内，如果在该范围内，则将对应的 yy 值置为0，表示该频率的信号被滤除；
    if m*(fs/n)<f1 || m*(fs/n)>f2 || m*(fs/n)<62 && m*(fs/n) > 60  %将奈奎斯特之后的频率也滤除点掉
        yy(m+1)=0;
    else
%         如果当前采样点对应的频率不在 8Hz 到 15Hz 范围内，则将 yy 的值保持为原始信号 y 的值。
        yy(m+1)=y(m+1);
    end
end %将频率为8Hz-15Hz的信号的幅值置0
yyi=abs(yy);

%%
% 滤波后
figure(3)
plot(x,yyi(1:n/2))
yi=ifft(yy);
figure(4)
subplot(212)
plot(t, real(yi), 'k') % 绘制滤波后的
xlabel('Time (\mu s)') % 设置 x 轴标签单位为微秒
ylabel('Amplitude')
title('filtered Signal')
