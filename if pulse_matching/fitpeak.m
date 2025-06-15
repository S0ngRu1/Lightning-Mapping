function fitted_peak_x = fitpeak(data,peak_index)
if peak_index+10 < 10240 && peak_index-10 > 0
    fit_range = (peak_index + (-10:10))';
elseif peak_index+6 < 10240 && peak_index-6 > 0
    fit_range = (peak_index + (-6:6))';
elseif peak_index+2 < 10240 && peak_index-2 > 0
    fit_range = (peak_index + (-2:2))';
else
    fitted_peak_x = peak_index;
    return;
end
fit_values = data(fit_range);
coefficients = polyfit(fit_range, fit_values, 2);
fit_indices_curve = linspace(min(fit_range), max(fit_range), 1000);
fit_values_curve = polyval(coefficients, fit_indices_curve);
% 绘制原始数据和拟合曲线
% figure;
% plot(1:length(data),data)
% plot(1:length(data),data, 'b', fit_indices_curve, fit_values_curve, 'r--');
% legend('原始数据', '拟合曲线');
% xlabel('y的索引');
% ylabel('y的值');
[~, max_index_fit] = max(fit_values_curve);
fitted_peak_x = fit_indices_curve(1,max_index_fit);
end
