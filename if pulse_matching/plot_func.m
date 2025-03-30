%���1
% logicalIndex = abs(result1.t123) < 1 & abs(result1.Rcorr) > 0.3;
logicalIndex = abs(result1.t123) < 1 & abs(result1.Rcorr) > 0.3 & result1.Start_loc < 523330212 & result1.Start_loc > 451518507;
filteredTable1 = result1(logicalIndex, :);
Start_loc = filteredTable1.Start_loc;
colorValues = (Start_loc - 3e8) / 2e8;

% 将方位角从0-360转换为-180到180
% filteredTable1.Azimuth = mod(filteredTable1.Azimuth + 180, 360) - 180;

figure;
scatter(filteredTable1.Azimuth, filteredTable1.Elevation, 1, colorValues, 'filled');
title('Azimuth vs Elevation');
xlabel('Azimuth');
% xlim([-180, 180]); % 修改x轴范围
% xticks(-180:40:180); % 修改x轴刻度
xlim([0, 360]);
xticks(0:40:360);
ylabel('Elevation');
ylim([-40, 100]);
yticks(-40:20:100);
colormap('hot');
colorbar;
caxis([0, 1.5]);
grid on;

%���2
figure;
% logicalIndex = abs(result2.t123) < 1 & abs(result2.Rcorr) > 0.3 &  result2.Start_loc > 500000000 ;
logicalIndex = abs(result2.t123) < 1 & abs(result2.Rcorr) > 0.5;
filteredTable2 = result2(logicalIndex, :);
index = 1:256:size(filteredTable2, 1);
Start_loc = filteredTable2.Start_loc;
colorValues = (Start_loc - 3e8) / 2e8;
figure;
scatter(filteredTable2.Azimuth,filteredTable2.Elevation, 1, colorValues, 'filled');
title('Azimuth vs Elevation');
xlabel('Azimuth');
xlim([0, 360]);
xticks(0:40:360);
ylabel('Elevation');
ylim([-40, 100]);
yticks(-40:20:100);
colormap('hsv');
colorbar;
caxis([0, 1.5]);
grid on;

%���3
logicalIndex = abs(result3.t123) > 0.001 & abs(result3.Rcorr) < 0.2 &  result3.Start_loc < 600000000 & result3.Start_loc > 400000000;
filteredTable3 = result3(logicalIndex, :);
Start_loc = filteredTable3.Start_loc;
colorValues = (Start_loc - 3e8) / 2e8;
figure;
scatter(filteredTable3.Azimuth,filteredTable3.Elevation, 1, colorValues, 'filled');
title('Azimuth vs Elevation');
xlabel('Azimuth');
xlim([0, 360]);
xticks(0:40:360);
ylabel('Elevation');
ylim([-40, 100]);
yticks(-40:20:100);
colormap('hsv');
colorbar;
caxis([0, 1.5]);
grid on;































% logicalIndex = abs(Untitled.VarName11) < 1 & abs(Untitled.VarName12) > 0.25 & Untitled.Paramers >500000000;
logicalIndex = abs(Untitled.VarName11) < 1 & abs(Untitled.VarName12) > 0.3 & Untitled.Paramers > 400000000 & Untitled.Paramers < 600000000 ;
filteredTable3 = Untitled(logicalIndex, :);
Start_loc = filteredTable3.Paramers;
colorValues = (Start_loc - 3e8) / 2e8;
figure;
scatter(filteredTable3.VarName5,filteredTable3.VarName6, 1, colorValues, 'filled');
title('Azimuth vs Elevation');
xlabel('Azimuth');
xlim([0, 360]);
xticks(0:40:360);
ylabel('Elevation');
ylim([-40, 100]);
yticks(-40:20:100);
colormap('hot');
colorbar;
caxis([0, 1.5]);
grid on;



% logicalIndex = abs(result1.t123) < 1 & abs(result1.Rcorr) > 0.3;
% filteredTable1 = result1(logicalIndex, :);
% index = 1:256:size(filteredTable1, 1);
% Start_loc = filteredTable1.Start_loc;
% colorValues = (Start_loc - 3e8) / 2e8;
% figure;
% scatter(filteredTable1.Azimuth,filteredTable1.Elevation, 1, colorValues, 'filled');
% title('Azimuth vs Elevation');
% xlabel('Azimuth');
% xlim([0, 360]);
% xticks(0:40:360);
% ylabel('Elevation');
% ylim([-40, 100]);
% yticks(-40:20:100);
% colormap('hot');
% colorbar;
% caxis([0, 1.5]);
% grid on;






% ������������
plot(processed_yld_signal-50, 'b');
hold on;
plot(processed_chj_signal, 'r');
legend('yld', 'chj');
xlabel('��������');
ylabel('��ֵ');
title('bp filter');






yld_signal = read_signal('../20240822165932.6610CH1.dat',6000,start_read_loc_yld);
chj_signal = read_signal('../2024 822 85933.651462CH1.dat',6000,start_read_loc_chj);
x = downsample(yld_signal,1);
y = downsample(chj_signal,1);
subplot(2,1,1);plot(x);title('yld');xlabel('��������');ylabel('��ֵ');
subplot(2,1,2);plot(y);title('chj');xlabel('��������');ylabel('��ֵ');




scatter(filteredTable1.cos_alpha_opt,filteredTable1.cos_beta_opt,2);
xlim([-1, 1]);
% ָ��y�᷶Χ�Ϳ̶ȱ��
ylim([-1, 1]);

scatter(Untitled.VarName5,Untitled.VarName6,1);
xlim([0, 360]);
xticks(0:40:360);
% ָ��y�᷶Χ�Ϳ̶ȱ��
ylim([0, 90]);
yticks(0:20:90);
scatter(result1.VarName10,result1.VarName11);
scatter(result1.cos,result1.cos1);
histogram(result.VarName11);
histogram(Untitled.VarName6);
subplot(3,1,1);plot(signal3);title('ch1');xlabel('��������');ylabel('��ֵ');
%���ú�������Ϊ32

subplot(3,1,2);plot(filtered_signal3);title('ch2');xlabel('��������');ylabel('��ֵ');

subplot(3,1,3);plot(filtered_signal3);title('ch3');xlabel('��������');ylabel('��ֵ');


plot(lag12_msw,R12_msw);
plot(ch1);
 %�����ϲ����Ա�ͼ
plot(ch3_gcc_new, 'b');
hold on;
plot(ch3_upsp(:,1), ch3_upsp(:,2), 'r--');
legend('ch3����', '�ϲ���ch3����');
xlabel('��������');
ylabel('��ֵ');
plot(lag12_gcc,R12_gcc);

 %�����ϲ����Ա�ͼ
plot(filtered_signal1, 'b');
hold on;
plot(ch1, 'r--');
legend('ch2����', '�ϲ���ch2����');
xlabel('��������');
ylabel('��ֵ');
axis auto


% ������������
plot(ch1, 'b');
hold on;
plot(IMF(:,3), 'r--');
legend('ch1����', '���˺��ch1����');
xlabel('��������');
ylabel('��ֵ');
title('���ݻ�������ߵ����ֵ��������ƽ�ƵĽ��');

% ������������
plot(ch1_gcc_new, 'b');
hold on;
plot(ch3_gcc_new, 'r--');
legend('ԭ�ź�', 'ƽ�ƺ�');
xlabel('��������');
ylabel('��ֵ');
title('���');



num_bins = 100; % ����ֱ��ͼ������
[counts, edges] = histcounts(matched_max_r_gccs, num_bins);

% ����ÿ����������ĵ㣬���ڻ���Ƶ�ʷֲ�ͼ
centers = edges(1:end-1) + diff(edges) / 2;

% ����Ƶ�ʷֲ�ͼ
figure;
bar(centers, counts);
xlabel('���ϵ����С');
ylabel('����');
title('ƥ�䵽�����ϵ���ֲ�ͼ');
grid on; % ��������

current_chj_start_loc = current_chj_read_loc + subsignal_starts(6) + floor(all_t_gccs(6));
chj_signal = read_signal('../2024 822 85933.651462CH1.dat', yld_signal_length, current_chj_start_loc);
plot(chj_signal)
figure;
plot(yld_signal)