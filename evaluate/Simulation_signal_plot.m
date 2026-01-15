%% Signal Processing Workflow Visualization for JGR Publication
% Updates:
% 1. Reduced to 2 subplots.
% 2. Shared Y-axis label for the tiled layout.
clear; clc; close all;

%% === 1. Signal Generation ===
fs = 200e6;                     
N = 4096;                       
win_len = 1024;                 
snr_val = 15;                   
true_delay_samples = 20.5;      
compress_scale = 0.4; % Scaling factor for raw signal visualization

% Filter Design: 20-80 MHz Bandpass
f_low = 20e6; f_high = 80e6;

% --- Generate Bipolar Pulse Cluster ---
rng(42); 
sig = zeros(1, N);
num_pulses = 40; 
range_start = floor(N/2) - 400;
range_end = floor(N/2) + 400;
pulse_locs = sort(randi([range_start, range_end], 1, num_pulses));

for p = 1:num_pulses
    width = randi([5, 12]);
    pulse_sign = sign(randn()); 
    amp = pulse_sign * (0.5 + 0.5*rand());
    p_start = max(1, pulse_locs(p)-width);
    p_end = min(N, pulse_locs(p)+width);
    L = p_end - p_start + 1;
    sig(p_start:p_end) = sig(p_start:p_end) + amp * hann(L)'; 
end
sig = sig / max(abs(sig));

% --- Apply Delay and Noise ---
f_axis = [0:ceil(N/2)-1, -floor(N/2):-1] * (fs/N);
phase_shift = exp(-1j * 2 * pi * f_axis * (true_delay_samples/fs));
SIG = fft(sig);
sig_delayed = ifft(SIG .* phase_shift, 'symmetric');

sig_power = mean(sig.^2);
noise_power = sig_power / (10^(snr_val/10));
noise_std = sqrt(noise_power);

% Generate Raw Noisy Signals (Original & Delayed)
r1 = sig + noise_std * randn(1, N);
r2 = sig_delayed + noise_std * randn(1, N);

% --- Filtering (Zero-Phase) and Windowing ---
r1_filt = filter_bp(r1, f_low, f_high, 5);
r2_filt = filter_bp(r2, f_low, f_high, 5);

center_idx = floor(N/2);
w_start = center_idx - floor(win_len/2);
idx = w_start : (w_start + win_len - 1);

ham_win = hamming(win_len)';
w1_final = r1_filt(idx) .* ham_win;

%% === 2. JGR Format Plotting ===
% Height adjusted for 2 panels
figure('Units', 'inches', 'Color', 'w'); 
t = tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

% Define Plot Range
plot_range = (w_start - 200) : (w_start + win_len + 200);
plot_range = plot_range(plot_range > 0 & plot_range <= N);

% JGR Style Settings
fontName = 'Arial'; 
fontSize = 11;
lineWidthMain = 1.2;
lineWidthThin = 0.8;

% === Panel (a): Original vs Delayed (Raw) ===
ax1 = nexttile; hold on; box on;
% Plot
plot(plot_range, r1(plot_range) * compress_scale, 'Color', [0.7 0.7 0.7], 'LineWidth', lineWidthThin, ...
    'DisplayName', 'Original Signal');
plot(plot_range, r2(plot_range) * compress_scale, 'Color', [0.85 0.325 0.098], 'LineWidth', lineWidthThin, ...
    'DisplayName', 'Delayed Signal');

% [修改点] 移除了单独的 ylabel
xlim([plot_range(1), plot_range(end)]);
ylim([-0.5, 0.5]); 
set(gca, 'XTickLabel', []); 
grid on; set(gca, 'GridLineStyle', ':');

legend('Location', 'northeast', 'Box', 'off', 'Interpreter', 'tex', 'FontName', fontName);
text(0.02, 0.9, '(a)', 'Units', 'normalized', 'FontName', fontName, 'FontSize', 12, 'FontWeight', 'bold');

% === Panel (b): Windowing Process (Formerly c) ===
ax2 = nexttile; hold on; box on;
% Filtered (Pre-window background)
plot(plot_range, r1_filt(plot_range), '--', 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 1, ...
    'DisplayName', 'Filtered (Pre-window)');
% Windowed (Final)
plot(idx, w1_final, 'Color', [0 0 0.8], 'LineWidth', 1.5, ...
    'DisplayName', 'Windowed Signal');
% Hamming Window visualization
plot(idx, ham_win * max(abs(w1_final))*1.2, 'Color', [0.4660 0.6740 0.1880], 'LineStyle', ':', 'LineWidth', 1.5, ...
    'DisplayName', 'Hamming Window');

xlabel('Time (samples)', 'FontName', fontName, 'FontSize', fontSize);

% [修改点] 移除了单独的 ylabel
xlim([plot_range(1), plot_range(end)]);
ylim([-0.1, 0.2]); 
grid on; set(gca, 'GridLineStyle', ':');

% Visual Boundary
xline(idx(1), 'k-', 'LineWidth', 0.5, 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
xline(idx(end), 'k-', 'LineWidth', 0.5, 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');

legend('Location', 'northeast', 'Box', 'off', 'Interpreter', 'tex', 'FontName', fontName);
text(0.02, 0.9, '(b)', 'Units', 'normalized', 'FontName', fontName, 'FontSize', 12, 'FontWeight', 'bold');

% === [关键修改] 添加共享 Y 轴标题 ===
ylabel(t, 'Amplitude', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'normal');

% Global Font Adjustment
set(findall(gcf,'-property','FontSize'),'FontSize', fontSize);
set(findall(gcf,'-property','FontName'),'FontName', fontName);
fprintf('Figure updated: Shared Y-label added.\n');

%% Helper Function
function filtered_signal = filter_bp(signal,f1,f2,order)
    Fs = 200e6;
    fn = Fs/2;
    Wn = [f1 f2]/fn;
    [b,a] = butter(order,Wn);
    filtered_signal = filtfilt(b,a,signal);
end