%�1�7�1�7�1�71
% logicalIndex = abs(result1.t123) <1 & abs(result1.Rcorr) > 0.3 &  result1.Start_loc < 4.8e8 & result1.Start_loc > 4.68e8;
logicalIndex = abs(result1.t123) <1 & abs(result1.Rcorr) > 0.3 &  result1.Start_loc < 4.705+34151156 & result1.Start_loc > 4.698e8+34151156;
filteredTable1 = result1(logicalIndex, :);
Start_loc = filteredTable1.Start_loc;
% colorValues = (Start_loc - 3e8) / 2e8;
colorValues = (Start_loc - min(Start_loc)) / (max(Start_loc) - min(Start_loc));
% �0�8�1�7�1�7 Azimuth �1�7�1�7���1�7�1�7 0-360 �1�7�1�7 -180-180
% filteredTable1.Azimuth = mod(filteredTable1.Azimuth - 180, 360) - 180;
figure;
scatter(filteredTable1.Azimuth,filteredTable1.Elevation, 2, colorValues, 'filled');
title('Azimuth vs Elevation');
xlabel('Azimuth');
xlim([0, 360]);
xticks(0:40:360);
ylabel('Elevation');
ylim([0, 100]);
yticks(0:10:90);
colormap('hot');
colorbar;
caxis([0, 1.5]);
grid on;
% 
% %�1�7�1�7�1�72
% figure;
%  logicalIndex = abs(result1.t123) <1 & abs(result1.Rcorr) > 0.4 &  result1.Start_loc < 4.5e8 & result1.Start_loc > 4e8;
logicalIndex = abs(result2.t123) < 1 & abs(result2.Rcorr) > 0.3;
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
colormap('jet');
colorbar;
caxis([0, 1.5]);
grid on;

% %�1�7�1�7�1�73
% logicalIndex = abs(result3.t123) > 0.001 & abs(result3.Rcorr) < 0.2 &  result3.Start_loc < 600000000 & result3.Start_loc > 400000000;
% filteredTable3 = result3(logicalIndex, :);
% Start_loc = filteredTable3.Start_loc;
% colorValues = (Start_loc - 3e8) / 2e8;
% figure;
% scatter(filteredTable3.Azimuth,filteredTable3.Elevation, 1, colorValues, 'filled');
% title('Azimuth vs Elevation');
% xlabel('Azimuth');
% xlim([0, 360]);
% xticks(0:40:360);
% ylabel('Elevation');
% ylim([-40, 100]);
% yticks(-40:20:100);
% colormap('hsv');
% colorbar;
% caxis([0, 1.5]);
% grid on;

signal1 = downsample(ch1,50);
signal2 = downsample(ch2,50);
signal3 = downsample(ch3,50);
plot_signal_spectrum(ch1(1:7000000));
plot_signal_spectrum(ch1(1e7:1.5e7));
plot_signal_spectrum(ch1(2.4e7:2.6e7));

data_original1=window_plus(N,data_original1);
data_original2=window_plus(N,data_original2);
Wc=2*80e6/200e6;
Wp=2*30e6/200e6;
data_filter1 = datafilter(data_original1,Wc,Wp);
data_filter2 = datafilter(data_original2,Wc,Wp);
data_noisefilter1 = datafilter(data_noise1,Wc,Wp);
data_noisefilter2 = datafilter(data_noise2,Wc,Wp);


Wc=2*80e6/200e6;
Wp=2*30e6/200e6;
N = 2e7;
r_loction = 256;
ch1 = read_signal('..\\20240822165932.6610CH1.dat',N,r_loction);
data_original1=window_plus(N,ch1);
data_filter1 = datafilter(data_original1,Wc,Wp);


signal_length = 1024;
r_loction_yld = 469401640;
r_loction_chj = 504138238;
ch_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction_yld-256);
ch_chj = read_signal('..\\2024 822 85933.651462CH1.dat',signal_length,r_loction_chj);
% filtered_yld = rfi_filter(ch_yld,1536);
% filtered_chj = rfi_filter(ch_chj,1536);
filtered_yld = filter_bp(ch_yld,20e6,80e6,5);
filtered_chj = filter_bp(ch_chj,20e6,80e6,5);
subplot(2,1,1);plot(filtered_yld);title('yld');xlabel('�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7');ylabel('�1�7�1�7�0�5');
subplot(2,1,2);plot(filtered_chj);title('chj');xlabel('�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7');ylabel('�1�7�1�7�0�5');

 filtered_chj_signal1 = waveletDenoiseAdaptive(ch2, level, wavelet);
 filtered_chj_signal2 = filter_bp(ch2,20e6,80e6,5);
x1 = downsample(ch1,50);
x2 = downsample(ch2,50);
x3 = downsample(ch3,50);
% x3 = downsample(ch3,50);
subplot(3,1,1);plot(ch1);title('ch1');xlabel('�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7');ylabel('�1�7�1�7�0�5');
subplot(3,1,2);plot(data_filter1);title('wf');xlabel('�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7');ylabel('�1�7�1�7�0�5');
subplot(3,1,3);plot(filter2);title('bp');xlabel('�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7');ylabel('�1�7�1�7�0�5');

%�1�7�1�7�1�7�1�725%�1�7�0�2�1�7


%�1�7�0�2�0�4�1�7�1�7�1�7�1�7�1�7�1�7�1�7�0�5�1�7�0�3�1�7�0�0
N = 5e7; % �1�7�0�2�1�7�1�7�1�7�1�7�1�7�1�7
subsignal_length = 4000;
M = ceil(N / subsignal_length); % �1�7�1�0�1�7�1�7�1�7
subsignal_start = 1:subsignal_length:length(filtered_signal1);
for subi = 1:numel(subsignal_start)
subsignal1 = filtered_signal1(subsignal_start(subi):subsignal_start(subi)+subsignal_length-1);
threshold =  mean(abs(subsignal1)) + 3*std(subsignal1);
end
% �1�7�1�7�1�7�1�7�0�7�1�7�ń1�7�1�7�1�7�0�5�1�7�0�4�1�7�1�7�1�7�1�7�1�7�0�0�1�7�1�7���1�7�1�7 thresholds �1�7�1�7
threshold = rand(1, M) * 5; % �0�5�1�7�1�7�1�7�1�7�0�5�1�7�1�7�0�6�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7�0�5�1�7�I
% �1�7�1�7�1�7�1�7�1�7�0�2�1�7
figure;
plot(filtered_signal1, 'b', 'LineWidth', 1.5); % �1�7�1�7�1�7�1�7�1�7�0�2�0�0�1�7�1�7�1�7�0�2�1�7�1�7
hold on;
% �1�7�1�7�1�7�1�7�0�7�1�7�ń1�7�1�7�1�7�0�5
for i = 1:M
    start_index = (i - 1) * subsignal_length + 1;
    end_index = min(i * subsignal_length, N); % �1�7�1�7�0�9�1�7�1�7�1�7�0�5�1�7�Ä1�7�1�7�1�7�1�7�0�2�0�5�1�7�1�7�1�7
    x = [start_index, end_index]; % �1�7�1�7�0�2�1�7�ń1�7 x �1�7�1�7��
    y = [threshold(i), threshold(i)]; % �1�7�1�7�0�2�1�7�ń1�7�1�7�1�7�0�5
    plot(x, y, 'r--', 'LineWidth', 1.5); % �1�7���1�7�0�2�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7�0�5
end
% �1�7�1�7�1�7�1�7�0�0�1�7�1�7�1�7�0�9�1�7�0�5
legend('filtered_signal1', 'Thresholds');
xlabel('Sample Index');
ylabel('Amplitude');
title('Signal and Thresholds Comparison');
grid on; % �1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7
hold off;





%3D�0�82D
x = filtered_S(:, 1);
y = filtered_S(:, 2);
z = filtered_S(:, 3);
% �1�7�1�7�0�1�1�7�1�7�1�7�1�7�1�7�1�7�0�8�1�7�1�7�0�2�1�7�1�7�1�7�1�7�1�7�1�7
[azimuth, elevation, ~] = cart2sph(x, y, z);
% �1�7�1�7�1�7�1�7�1�7�1�7�0�8�1�7�1�7�0�2�1�7�0�8�1�7
azimuth = azimuth * 180 / pi;
elevation = elevation * 180 / pi;
% �1�7�1�7�1�7�1�7�1�7�1�7�˄1�7�0�9�1�7���0�2 0-360 �1�7�1�7
azimuth = mod(azimuth, 360);
azimuth(azimuth < 0) = azimuth(azimuth < 0) + 360;
% �1�7�1�7�1�7�1�7�1�7�1�7�1�7�0�9�1�7���0�2 0-90 �1�7�1�7
elevation = abs(elevation); % �0�2�1�7�1�7�1�7�1�7�1�7�1�7�0�2�1�7�0�0�1�7�0�5
elevation(elevation > 90) = 90; % �1�7�1�7�1�7�1�7�1�7�1�7�1�7�0�5�0�2 90 �1�7�1�7
% �1�7�1�7�1�7�0�3�1�7�˄1�7�0�1�1�7�1�7�1�7�1�7�0�7�1�7�0�3�1�7�1�7�0�0
figure;
scatter(azimuth, elevation, 1, 'filled');
xlabel('�1�7�1�7�˄1�7�1�7 (�1�7�1�7)');
ylabel('�1�7�1�7�1�7�1�7 (�1�7�1�7)');
title('�1�7�1�7�˄1�7�0�1�1�7�1�7�1�7�1�7�0�9�0�2�1�7');
grid on;
axis([0 360 0 90]); % �1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7�5�4��
colorbar; % �1�7�1�7�1�7�1�7�1�7�1�7�0�2�1�7�1�7�1�7�1�7�1�7�1�7�0�5�1�7�1�7




subplot(3,1,1);plot(ch1);title('ch1');xlabel('�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7');ylabel('�1�7�1�7�0�5');
subplot(3,1,2);plot(ch2);title('ch2');xlabel('�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7');ylabel('�1�7�1�7�0�5');
subplot(3,1,3);plot(ch3);title('ch3');xlabel('�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7');ylabel('�1�7�1�7�0�5');







fs = 200e6;
signal_length = 71916308;
r_loction = 451501808;
ch1_yld = read_signal('..\\20240822165932.6610CH1.dat',signal_length,r_loction);
ch1_chj = read_signal('..\\2024 822 85933.651462CH1.dat',signal_length,r_loction + 3.3e7);
ch2 = read_signal('..\\2024 822 85933.651462CH2.dat',signal_length,r_loction);
ch3 = read_signal('..\\2024 822 85933.651462CH3.dat',signal_length,r_loction+165/5);

% �1�7�1�7�0�0�1�7�0�0�1�7�0�2�1�7�1�7�1�7
[b, a] = butter(4, [20e6, 80e6] / (fs / 2), 'bandpass'); 
% �1�7�1�7�1�7�0�2�0�5�1�7�1�7�Մ1�7�0�0�1�7�0�2�1�7
filtered_signal1 = filter(b, a, ch1);
filtered_signal2 = filter(b, a, ch2);
filtered_signal3 = filter(b, a, ch3);

%bp
filtered_signal1_yld = filter_bp(ch1_yld, 20e6 ,80e6 ,5);
filtered_signal1_chj = filter_bp(ch1_chj, 20e6 ,80e6 ,5);
filtered_signal3 = filter_bp(ch3, 20e6 ,80e6 ,5);








subplot(3,1,1);plot(filtered_signal1);title('ch1');xlabel('�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7');ylabel('�1�7�1�7�0�5');
subplot(3,1,2);plot(filtered_signal2);title('ch2');xlabel('�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7');ylabel('�1�7�1�7�0�5');
subplot(3,1,3);plot(filtered_signal3);title('ch3');xlabel('�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7');ylabel('�1�7�1�7�0�5');


plot(ch3,'b');
hold on;
plot(filtered_signal3 +50,'r');
legend('�0�9�1�7�0�2�1�7','�1�7�0�2�1�7�1�7�1�7');
xlabel('�1�7�1�7�1�7�1�7�1�7�1�7�1�7�1�7');
ylabel('�1�7�1�7�0�5');


[R1_x, R1_y, R1_z] = sph2cart(deg2rad(181.588691),deg2rad(49.691292),1);
% �1�7�1�7�1�7�1�7�1�7�1�7�1�7�0�9�1�7�1�7�1�7�˄1�7�1�7�1�7�0�2�1�7
theta = asin(R1_z) * (180/pi); % �0�0�0�5�1�7�1�7���0�2 [-90, 90] �1�7�1�7
% �1�7�1�7�1�7�1�2�˄1�7�0�9�1�7�1�7�1�7�˄1�7�1�7�1�7�0�2�1�7
phi = atan2(R1_y, R1_x) * (180/pi); % �0�0�0�5�1�7�1�7���0�2 [-180, 180] �1�7�1�7
% �1�7�1�7�1�7�1�7�1�7�0�8�1�7�1�7�1�7�1�7�˄1�7�1�7�0�8�1�7�1�7�0�2 [0, 360) �1�7�1�7��
if phi < 0
    phi = phi + 360;
end
% �1�7�1�7�1�7�1�7�1�7�1�7
fprintf('�1�7�1�7�˄1�7�1�7 (Azimuth): %.2f �1�7�1�7\n', phi);
fprintf('�1�7�1�7�1�7�1�7 (Elevation): %.2f �1�7�1�7\n', theta);



