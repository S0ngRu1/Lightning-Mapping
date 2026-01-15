%% Lightning Location Algorithm Assessment: JGR Publication Quality
% Sampling Rate: 200 MHz
% Updates: Arial fonts, English labels, distinct markers, panel labels (a)/(b)

clear; clc; close all;

%% === 1. Global Simulation Parameters ===
fs = 200e6;                 
T_sim = 50e-6;              
N_total = round(T_sim * fs); 
true_delay_samples = 2.4;   
num_trials = 300;           
SNR_dB_range = -10:5:15; % Adjusted to match the uploaded image range roughly
f1 = 20e6; f2 = 80e6;
[b_filt, a_filt] = butter(4, [f1, f2]/(fs/2)); 

fprintf('Simulation Started (Parfor enabled)...\n');

%% ========================================================================
%  Experiment 1: Upsampling Strategy Comparison (Figure 1 data)
% ========================================================================
fprintf('Running Exp 1: Upsampling Strategies...\n');
fixed_win = 1024;
fixed_sig_len = 800; 
test_upsamples = [1, 5, 10, 30, 50]; 
corr_interp_compare = 8; 

% 1. Strategy A: Raw Signal Upsampling
res_sig_up = zeros(length(test_upsamples), length(SNR_dB_range));
for i = 1:length(test_upsamples)
    curr_up = test_upsamples(i);
    fprintf('  -> [Strategy A] Raw Signal Upsampling K=%d\n', curr_up);
    res_sig_up(i, :) = run_core_sim(SNR_dB_range, num_trials, fixed_win, ...
                                    curr_up, 0, ...
                                    fixed_sig_len, fs, N_total, true_delay_samples, b_filt, a_filt);
end

% 2. Strategy B: Correlation Interpolation (Control)
fprintf('  -> [Strategy B] Correlation Interp M=%d (Control)\n', corr_interp_compare);
res_corr_up = run_core_sim(SNR_dB_range, num_trials, fixed_win, ...
                           1, corr_interp_compare, ...
                           fixed_sig_len, fs, N_total, true_delay_samples, b_filt, a_filt);

%% ========================================================================
%  Experiment 2: Adaptive Windowing (Figure 2 data)
% ========================================================================
fprintf('Running Exp 2: Adaptive Windowing...\n');
scenarios = [400, 900, 1900, 3900];
scenario_names = {'High Density', 'Medium Density', 'Low Density', 'Sparse Signal'};
test_wins = [512, 1024, 2048, 4096];
fixed_upsample = 10;
res_exp2 = cell(1, 4); 

for s = 1:length(scenarios)
    curr_len = scenarios(s);
    fprintf('  -> Scenario %d: Signal Len %d\n', s, curr_len);
    scene_res = zeros(length(test_wins), length(SNR_dB_range));
    for w = 1:length(test_wins)
        curr_win = test_wins(w);
        scene_res(w, :) = run_core_sim(SNR_dB_range, num_trials, curr_win, ...
                                       fixed_upsample, 0, ...
                                       curr_len, fs, N_total, true_delay_samples, b_filt, a_filt);
    end
    res_exp2{s} = scene_res;
end

%% ========================================================================
%  JGR Formatting & Plotting
% ========================================================================
% General JGR Settings
fontName = 'Arial';
fontSize = 10;
lineWidthMain = 1.5;
lineWidthThin = 1.0;
markerSize = 6;

% --- Figure 1: Upsampling Strategy Comparison ---
figure('Units', 'inches', 'Position', [1, 1, 7, 5], 'Color', 'w'); % Standard 1-column or 2-column width
ax1 = gca;
hold on; box on; grid on;

% Plot Strategy A (Colored Lines)
% Use distinct colors suitable for publication
colors = [0 0.4470 0.7410;  % Blue
          0.3010 0.7450 0.9330; % Cyan
          0.4660 0.6740 0.1880; % Green
          0.9290 0.6940 0.1250; % Yellow
          0.8500 0.3250 0.0980]; % Orange

legend_str = {};
for i = 1:length(test_upsamples)
    plot(SNR_dB_range, res_sig_up(i, :), 'o-', 'LineWidth', lineWidthMain, ...
         'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', markerSize);
    legend_str{end+1} = sprintf('Raw Upsample (K=%d)', test_upsamples(i));
end

% Plot Strategy B (Black Dashed)
plot(SNR_dB_range, res_corr_up, 's--', 'Color', 'k', 'LineWidth', lineWidthMain, ...
     'MarkerFaceColor', 'k', 'MarkerSize', markerSize);
legend_str{end+1} = sprintf('Corr. Interp. (M=%d)', corr_interp_compare);

% Labels & Scale
xlabel('SNR (dB)', 'FontName', fontName, 'FontSize', fontSize);
ylabel('Time Delay STD (s)', 'FontName', fontName, 'FontSize', fontSize);
set(gca, 'YScale', 'log', 'FontName', fontName, 'FontSize', fontSize);
ylim([1e-11, 1e-6]); 
grid on; set(gca, 'GridLineStyle', ':');

% Legend
legend(legend_str, 'Location', 'northeast', 'Box', 'off', 'FontName', fontName, 'FontSize', 9);

% --- Figure 2: Adaptive Windowing ---
figure('Units', 'inches', 'Position', [1, 1, 7, 6], 'Color', 'w');
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for s = 1:4
    nexttile;
    hold on; box on; grid on;
    data = res_exp2{s};
    
    % Plot lines
    styles = {'-', '--', '-.', ':'}; % Distinct line styles
    for w = 1:4
        plot(SNR_dB_range, data(w, :), 'LineStyle', styles{w}, ...
            'LineWidth', 1.2, 'Marker', '.', 'MarkerSize', 12);
    end
    
    % Panel formatting
    title(scenario_names{s}, 'FontName', fontName, 'FontSize', 10, 'FontWeight', 'bold');
    if s > 2, xlabel('SNR (dB)', 'FontName', fontName, 'FontSize', fontSize); end
    if mod(s,2)==1, ylabel('STD (s)', 'FontName', fontName, 'FontSize', fontSize); end
    
    set(gca, 'YScale', 'log', 'FontName', fontName, 'FontSize', fontSize);
    ylim([1e-11, 1e-6]);
    grid on; set(gca, 'GridLineStyle', ':');
    
    % Add Panel Labels (b1, b2, etc. or just a,b,c,d)
    labels = {'(a)', '(b)', '(c)', '(d)'};
    text(0.05, 0.9, labels{s}, 'Units', 'normalized', 'FontName', fontName, 'FontSize', 11, 'FontWeight', 'bold');
    
    if s == 1
         legend({'Win=512', 'Win=1024', 'Win=2048', 'Win=4096'}, ...
               'Location', 'southwest', 'Box', 'off', 'FontSize', 8); 
    end
end

fprintf('Plots generated successfully.\n');

%% ========================================================================
%  Core Simulation Function (Unchanged logic, just keeping it functional)
% ========================================================================
function results = run_core_sim(SNRs, trials, win_len, sig_upsample, corr_upsample, sig_len, fs, N, true_delay, b, a)
    results = zeros(1, length(SNRs));
    f_axis = [0:ceil(N/2)-1, -floor(N/2):-1] * (fs/N);
    phase_shift = exp(-1j * 2 * pi * f_axis * (true_delay/fs));
    
    for i = 1:length(SNRs)
        snr_val = SNRs(i);
        errors = zeros(1, trials);
        parfor k = 1:trials
            raw_sig = randn(1, N);
            mask = zeros(1, N);
            center = floor(N/2);
            start_p = max(1, center - floor(sig_len/2));
            end_p = min(N, start_p + sig_len - 1);
            mask(start_p:end_p) = 1;
            sig = raw_sig .* mask;
            
            SIG = fft(sig);
            sig_delayed = ifft(SIG .* phase_shift, 'symmetric');
            
            current_sig_power = sum(sig.^2) / sig_len; 
            if current_sig_power == 0, current_sig_power = 1; end
            target_sig_power = 1 * 10^(snr_val/10); 
            scale = sqrt(target_sig_power / current_sig_power);
            
            r1 = sig * scale + randn(1, N);
            r2 = sig_delayed * scale + randn(1, N);
            
            r1 = filter(b, a, r1);
            r2 = filter(b, a, r2);
            w_start = floor(N/2) - floor(win_len/2);
            idx = w_start : (w_start + win_len - 1);
            idx = idx(idx>0 & idx<=N);
            w1 = r1(idx);
            w2 = r2(idx);
            
            if sig_upsample > 1
                x_old = 1:length(w1);
                x_new = linspace(1, length(w1), length(w1)*sig_upsample);
                w1_proc = interp1(x_old, w1, x_new, 'spline');
                w2_proc = interp1(x_old, w2, x_new, 'spline');
                lag_unit_scale = 1 / sig_upsample; 
            else
                w1_proc = w1; w2_proc = w2;
                lag_unit_scale = 1;
            end
            
            [cc, lags] = xcorr(w1_proc, w2_proc);
            
            if corr_upsample > 0
                cc_len = length(cc);
                c_old = 1:cc_len;
                c_new = linspace(1, cc_len, cc_len * corr_upsample);
                cc_final = interp1(c_old, cc, c_new, 'spline');
                lags_final = linspace(lags(1), lags(end), length(cc_final));
                lags_final = lags_final * lag_unit_scale;
            else
                cc_final = cc;
                lags_final = lags * lag_unit_scale;
            end
            
            [~, max_idx] = max(cc_final);
            
            if max_idx > 1 && max_idx < length(cc_final)
                y1 = cc_final(max_idx-1); 
                y2 = cc_final(max_idx); 
                y3 = cc_final(max_idx+1);
                delta = 0.5 * (y1 - y3) / (y1 - 2*y2 + y3);
            else
                delta = 0;
            end
            
            grid_resolution = lags_final(2) - lags_final(1);
            est_delay_pts = lags_final(max_idx) + delta * grid_resolution;
            est_delay_sec = est_delay_pts / fs;
            errors(k) = abs(est_delay_sec - (true_delay/fs));
        end
        results(i) = std(errors);
    end
end