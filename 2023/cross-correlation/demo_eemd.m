%% 1.生成仿真信号
fs = 400;  %采样频率
t = 0:1/fs:0.75; %时间轴
x = sin(2*pi*4*t); %低频正弦信号
y = 0.5*sin(2*pi*120*t); %高频正弦信号
for i = 1:length(t) %将高频信号处理成间断性
    if mod(t(i),0.25)>0.11&&mod(t(i),0.25)<0.12
    else
        y(i) = 0;
    end
end
sig = x+y; %信号叠加
figure('color','white')
plot(t,sig,'k') %绘制原始信号
%% 2.EEMD分解图
Nstd = 0.2; %Nstd为附加噪声标准差与Y标准差之比
NE = 100;   %NE为对信号的平均次数
imf = pEEMD(sig,t,Nstd,NE);
% 画信号EEMD分解图
% 输入：
% y为待分解信号
% FsOrT为采样频率或采样时间向量，如果为采样频率，该变量输入单个值；如果为时间向量，该变量为与y相同长度的一维向量。如果未知采样频率，可设置为1
% Nstd为附加噪声标准差与Y标准差之比
% NE为对信号的平均次数
% 输出：
% imf为经eemd分解后的各imf分量值（公开版代码无法输出imf）
%% 3.EEMD分解及频谱图
imf = pEEMDandFFT(sig,fs,Nstd,NE);% 画信号EEMD分解与各IMF分量频谱对照图
% function imf = pEEMDandFFT(y,FsOrT,Nstd,NE)
% 输入：
% y为待分解信号
% FsOrT为采样频率或采样时间向量，如果为采样频率，该变量输入单个值；如果为时间向量，该变量为与y相同长度的一维向量
% Nstd为附加噪声标准差与Y标准差之比
% NE为对信号的平均次数
% 输出：
% imf为经eemd分解后的各imf分量值（公开版代码无法输出imf）


%% 4. 根据频率的直接筛选方法
indices = [4 , 5 , 6,7,8];  % 假设我们想要保留的是前两个IMF，需要注意res等同于最后一个IMF，在此例子中，res即IMF5，如果要保留第1、2个IMF和res，可以写成：indices = [1, 2 ,5];
filtered_imfs = imf(indices, :);
filtered_signal1 = sum(filtered_imfs, 1);
ori = sig;  %无噪声信号
fil = filtered_signal1;  %滤波后信号
figure('color','w')
subplot(311);plot(x,'k');title('原始信号（未加入噪声）')
subplot(312);plot(sig,'k');title('原始信号（加入噪声）')
subplot(313);plot(fil,'k');title('滤波后信号')