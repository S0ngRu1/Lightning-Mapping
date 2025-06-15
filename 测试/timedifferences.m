% 将时间数据转换为纳秒
timeInNanoseconds = timestamps(:,1) + ...          % 纳秒
                    timestamps(:,2) * 1e3 + ...     % 微秒
                    timestamps(:,3) * 1e6 + ...     % 毫秒
                    timestamps(:,4) * 1e9 + ...     % 秒
                    timestamps(:,5) * 60 * 1e9 + ...% 分
                    timestamps(:,6) * 3600 * 1e9;   % 时

% 计算相邻时间的差值
timeDifferences = diff(timeInNanoseconds);
timeDifferences = timeDifferences-1e9;
% 绘制时间差值图
figure;
plot(abs(timeDifferences), 'LineWidth', 1.5);
xlabel('样本索引');
ylabel('时间差 (纳秒)');
title('相邻时间差值图');
grid on;