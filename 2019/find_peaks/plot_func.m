
%计算t123
result1.t123 = result1.t12 + result1.t23-result1.t13;
%将结果中的相关系数设置为大于0.进行筛选
% logicalIndex = abs(result1.corr12) > 0.3 & abs(result1.corr13) > 0.3 & abs(result1.corr23) > 0.3 & abs(result1.t123) < 1;
logicalIndex = abs(result1.t123) < 1 & abs(result1.Rcorr) > 0.3;
filteredTable1 = result1(logicalIndex, :);
% %画图
scatter(filteredTable1.Azimuth,filteredTable1.Elevation,1);
% scatter(result1.Azimuth,result1.Elevation,1);
xlim([0, 360]);
xticks(0:40:360);
% 指定y轴范围和刻度标记
ylim([-40, 100]);
yticks(-40:20:100);

figure;
logicalIndex = abs(result2.t123) < 1 & abs(result2.Rcorr) > 0.3;
filteredTable2 = result2(logicalIndex, :);
% %画图
scatter(filteredTable2.Azimuth,filteredTable2.Elevation,1);
% scatter(result1.Azimuth,result1.Elevation,1);
xlim([0, 360]);
xticks(0:40:360);
% 指定y轴范围和刻度标记
ylim([-40, 100]);
yticks(-40:20:100);


scatter(filteredTable1.cos,filteredTable1.cos1,2);
xlim([-1, 1]);
% 指定y轴范围和刻度标记
ylim([-1, 1]);

figure;
logicalIndex = abs(Untitled.VarName11) < 0.001 & abs(Untitled.VarName12) > 0.3;
filteredTable2 = Untitled(logicalIndex, :);
scatter(filteredTable2.VarName5,filteredTable2.VarName6,1);
xlim([0, 360]);
xticks(0:40:360);
% 指定y轴范围和刻度标记
ylim([-40, 100]);
yticks(-40:20:100);

scatter(result1.VarName10,result1.VarName11);
scatter(result1.cos,result1.cos1);
histogram(result.VarName11);
histogram(Untitled.VarName6);

subplot(3,1,1);plot(ch1_gcc);title('ch1');xlabel('采样点数');ylabel('幅值');
%设置横坐标间隔为32
xticks(0:1024:3072);
subplot(3,1,2);plot(ch2_gcc);title('ch2');xlabel('采样点数');ylabel('幅值');
xticks(0:1024:3072);
subplot(3,1,3);plot(ch3_gcc);title('ch3');xlabel('采样点数');ylabel('幅值');
xticks(0:1024:3072);

plot(lag12_msw,R12_msw);
plot(ch1);
 %绘制上采样对比图
plot(ch3_gcc_new, 'b');
hold on;
plot(ch3_upsp(:,1), ch3_upsp(:,2), 'r--');
legend('ch3波形', '上采样ch3波形');
xlabel('采样点数');
ylabel('幅值');
plot(lag12_gcc,R12_gcc);

 %绘制上采样对比图
plot(filtered_signal1, 'b');
hold on;
plot(ch1, 'r--');
legend('ch2波形', '上采样ch2波形');
xlabel('采样点数');
ylabel('幅值');
axis auto

% 绘制两个波形
plot(ch1, 'b');
hold on;
plot(IMF(:,3), 'r--');
legend('ch1波形', '过滤后的ch1波形');
xlabel('采样点数');
ylabel('幅值');
title('根据互相关曲线的最大值进行数据平移的结果');

% 绘制两个波形
plot(ch1_gcc, 'b');
hold on;
plot(ch3_new(idx-98:idx+91), 'r--');
legend('ch1波形', '平移后ch3波形');
xlabel('采样点数');
ylabel('幅值');
title('根据互相关曲线的最大值进行数据平移的结果');