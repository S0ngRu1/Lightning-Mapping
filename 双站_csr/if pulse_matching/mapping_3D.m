clear
close all; % Close existing figures

%% Step 0: Configuration & Parameters ====================================
fprintf('Step 0: Configuring parameters...\n');

% --- File Paths & Data Loading ---
% Define paths to your data files
CHJ_FILTERED_CH1_PATH = '../chj_filtered_ch1_3-6.mat';
CHJ_FILTERED_CH2_PATH = '../chj_filtered_ch2_3-6.mat';
CHJ_FILTERED_CH3_PATH = '../chj_filtered_ch3_3-6.mat';
YLD_FILTERED_CH1_PATH = '../yld_filtered_ch1_3-6.mat';
YLD_RESULT_PATH = 'result_yld_th1_4_5e8-5_5e8.txt'; % Pre-calculated YLD event info

% --- Station & Geometry ---
YLD_SIT = [0, 0, 0];         % YLD station coordinates (meters) - Reference station
CHJ_SIT = [1991, -7841, 0];   % CHJ station coordinates (meters) - Relative to YLD? Ensure consistency.
P_BASELINE = CHJ_SIT - YLD_SIT; % Baseline vector from YLD to CHJ (meters)
% DIST_SITES = norm(P_BASELINE); % Calculated distance: sqrt(1991^2 + (-7841)^2) = 8090.0 m
% fprintf('Baseline distance: %.1f meters\n', DIST_SITES);

% --- Physical Constants ---
C_LIGHT = 299792458; % Speed of light in m/s - CORRECTED VALUE

% --- Processing Parameters ---
CHJ_SIGNAL_LENGTH = 1024;    % Length of window for fine correlation (samples)
MATCH_SIGNAL_LENGTH = 6000;    % Half-length of segment extracted around coarse match for refinement search (samples)
                               % Total refinement search length = 2 * MATCH_SIGNAL_LENGTH + 1

% --- Offsets & Timing ---
% Explanation of Offsets:
% CHJ_FILE_OFFSET: Seems to be a constant time/index difference between the start
%                  of the CHJ recording and the YLD recording reference time.
% YLD_OFFSET_3: An offset used when selecting the YLD segment for refinement correlation.
%               Needs clear justification - is it aligning features based on coarse match?
%               3e8 seems related to original signal read location in commented code.
%               If filtered signals are aligned, this might not be needed or should be calculated differently.
CHJ_FILE_OFFSET = 34371950;        % Offset added to CHJ raw data read (samples) - From commented Step0
SUB_FILTER_SIGNAL_LENGTH = 60000;  % Parameter from original code (used in YLD_OFFSET_3 calc)
YLD_OFFSET_3 = 3e8 + SUB_FILTER_SIGNAL_LENGTH / 4; % Offset for YLD refinement segment (samples) - NEEDS VERIFICATION/JUSTIFICATION

TGCC_TIME_MULTIPLIER = 5;        % Multiplier for t_gcc -> time (ns/sample). Assumes 5 ns sample time (200 Msps). VERIFY YOUR SAMPLE RATE!
TIME_VALIDATION_TOLERANCE_W = 30000; % Max allowed difference between geometric and measured TDOA (ns). 30000 ns = 30 us.

% --- Thresholds ---
YLD_PREFILTER_RCORR_THR = 0.3;   % Min Rcorr for initial YLD event selection
YLD_PREFILTER_T123_THR = 1;      % Max t123 for initial YLD event selection (What is t123?)
REFINE_CORR_THR = 0.15;          % Min Rcorr for accepting refined correlation peak

% --- YLD Event Selection ---
START_READ_LOC_YLD = 451518508; % Start index for reading YLD results
END_READ_LOC_YLD = 523330211;   % End index for reading YLD results

%% Step 1: Load Data ======================================================
fprintf('Step 1: Loading filtered data and YLD results...\n');

try
    load(CHJ_FILTERED_CH1_PATH, 'filtered_chj_signal1');
    load(CHJ_FILTERED_CH2_PATH, 'filtered_chj_signal2');
    load(CHJ_FILTERED_CH3_PATH, 'filtered_chj_signal3');
    load(YLD_FILTERED_CH1_PATH, 'filtered_yld_signal1');
    fprintf('Filtered data loaded successfully.\n');
catch ME
    fprintf('Error loading filtered data MAT files!\n');
    rethrow(ME);
end

% Read YLD pre-calculated results (location, azimuth, elevation, etc.)
try
    [yld_start_loc, yld_azimuth, yld_elevation, yld_Rcorr, yld_t123] = read_result(YLD_RESULT_PATH, START_READ_LOC_YLD, END_READ_LOC_YLD);
    fprintf('YLD results loaded successfully. %d events found.\n', numel(yld_start_loc));
catch ME
    fprintf('Error reading YLD result file: %s\n', YLD_RESULT_PATH);
    rethrow(ME);
end

% --- Data Length Checks ---
len_ch1 = length(filtered_chj_signal1);
len_ch2 = length(filtered_chj_signal2);
len_ch3 = length(filtered_chj_signal3);
len_yld1 = length(filtered_yld_signal1);
if ~(len_ch1 == len_ch2 && len_ch1 == len_ch3)
    warning('CHJ filtered signal lengths differ! Check synchronization/filtering.');
    % Potentially add error handling or use min length
end
fprintf('Data lengths: YLD1=%d, CHJ1=%d, CHJ2=%d, CHJ3=%d\n', len_yld1, len_ch1, len_ch2, len_ch3);


%% Step 2: Process Events (Matching, Localization, Validation) ===========
fprintf('Step 2: Starting event processing loop...\n');

% --- Initialization ---
S_results = []; % Stores final 3D locations for valid events
match_results = struct('yld_start_loc', {}, 'chj_loc_refined_start', {}, 'r_gccs_refined', {}, 'S_calc', {}, 'time_residual_ns', {}); % Detailed results
dltas_ns = []; % Stores time validation residuals (ns)

first_successful_coarse_chj_loc = 0;  % Track the first successful coarse match index
last_successful_coarse_chj_loc = -Inf; % Track the *previous* successful coarse match index

h = waitbar(0, 'Processing Events...');
num_events = numel(yld_start_loc);

for i = 1:num_events
    waitbar(i/num_events, h, sprintf('Processing Event %d/%d (YLD loc %d)', i, num_events, yld_start_loc(i)));

    % --- Prefiltering based on yld data ---
    if yld_Rcorr(i) < YLD_PREFILTER_RCORR_THR || yld_t123(i) > YLD_PREFILTER_T123_THR
        fprintf('Event %d: Skipping yld index %d due to low Rcorr (%.2f)\n', i, yld_start_loc(i), yld_Rcorr(i));
        continue
    end

    current_coarse_chj_loc = []; % Initialize coarse location for this iteration
    r_gccs_coarse = 0;

    % --- A. Coarse Matching Step (Find corresponding CHJ segment) ---
    if first_successful_coarse_chj_loc == 0
        % First valid event: Perform full coarse search
        search_hint_loc = 0; % No hint, search widely
        [temp_loc, r_gccs_coarse] = get_match_single_yld_chj_find_peak(filtered_chj_signal1, filtered_yld_signal1, yld_start_loc(i), search_hint_loc);
        if ~isempty(temp_loc)
            current_coarse_chj_loc = temp_loc;
            first_successful_coarse_chj_loc = current_coarse_chj_loc; % Store for future hints
            last_successful_coarse_chj_loc = current_coarse_chj_loc;
            fprintf('Event %d: First coarse match found at CHJ index %d (R=%.3f)\n', i, current_coarse_chj_loc, r_gccs_coarse);
        end
    else
        % Subsequent events: Use previous match(es) to guide search
        % Robust Option: Use *last* successful match as reference
        if i > 1 && ~isinf(last_successful_coarse_chj_loc) && i <= num_events && (i-1) > 0 && (i-1) <= num_events && yld_start_loc(i-1) > 0
             % Estimate based on previous step's relative offset in yld and chj
             search_hint_loc = last_successful_coarse_chj_loc + (yld_start_loc(i) - yld_start_loc(i-1));
             % Ensure search starts at least after the last match
             search_hint_loc = max(search_hint_loc, last_successful_coarse_chj_loc + 1);
             fprintf('Event %d: Using previous match hint. Search near CHJ index %d\n', i, search_hint_loc);
        else
            % Fallback: Use the first match as reference (less robust to drift)
             search_hint_loc = first_successful_coarse_chj_loc + (yld_start_loc(i) - yld_start_loc(1));
             search_hint_loc = max(search_hint_loc, last_successful_coarse_chj_loc + 1); % Ensure starts after last
             fprintf('Event %d: Using first match hint. Search near CHJ index %d\n', i, search_hint_loc);
        end

        [temp_loc, r_gccs_coarse] = get_match_single_yld_chj_find_peak(filtered_chj_signal1, filtered_yld_signal1, yld_start_loc(i), search_hint_loc);

        % **** STRICT NON-CROSSING ENFORCEMENT ****
        if ~isempty(temp_loc)
            if temp_loc > last_successful_coarse_chj_loc
                current_coarse_chj_loc = temp_loc;
                last_successful_coarse_chj_loc = current_coarse_chj_loc; % Update last successful loc
                fprintf('Event %d: Sequential coarse match found at CHJ index %d > %d (R=%.3f)\n', i, current_coarse_chj_loc, last_successful_coarse_chj_loc, r_gccs_coarse);
            else
                fprintf('Event %d: Warning! Coarse match %d <= previous %d. Skipping event (Non-sequential).\n', i, temp_loc, last_successful_coarse_chj_loc);
%                 current_coarse_chj_loc = []; % Discard non-sequential match
            end
        end
    end

    % --- Check if coarse match was successful and sequential ---
    if isempty(current_coarse_chj_loc)
        fprintf('Event %d: Failed to find a valid sequential coarse match for yld index %d\n', i, yld_start_loc(i));
        continue % Skip to the next yld segment
    end

     % --- Optional: Check coarse correlation quality ---
    % if r_gccs_coarse < SOME_THRESHOLD
    %     fprintf('Event %d: Skipping yld index %d due to low coarse correlation (%.2f)\n', i, yld_start_loc(i), r_gccs_coarse);
    %     continue;
    % end

    % --- B. Extract CHJ Signals for Refinement ---
    coarse_chj_loc = current_coarse_chj_loc; % Validated coarse location
    refine_search_start = max(1, coarse_chj_loc - MATCH_SIGNAL_LENGTH);
    refine_search_end = min(len_ch1, coarse_chj_loc + MATCH_SIGNAL_LENGTH); % Use actual CHJ length

    % Check if refinement window is valid and long enough
    if refine_search_start >= refine_search_end || (refine_search_end - refine_search_start + 1) < CHJ_SIGNAL_LENGTH
        fprintf('Event %d: Warning! Invalid refinement extraction window [%d, %d]. Skipping event.\n', i, refine_search_start, refine_search_end);
        continue;
    end

    chj_match_signal1 = filtered_chj_signal1(refine_search_start:refine_search_end);
    chj_match_signal2 = filtered_chj_signal2(refine_search_start:refine_search_end);
    chj_match_signal3 = filtered_chj_signal3(refine_search_start:refine_search_end);

    % --- C. Sliding Window Refinement & 3D Localization ---
    sub_S_results = []; % Store potential source locations for this event i
    sub_R_gccs = [];    % Store corresponding correlations
    sub_refined_indices = []; % Store start index of refined window
    sub_t_gcc_samples = []; % Store t_gcc from refinement

    % Define YLD segment for refinement correlation
    % ** CRITICAL: Verify YLD_OFFSET_3 logic. Is it necessary if signals are aligned? **
    % Assuming it's needed for now:
    yld_refine_start = yld_start_loc(i) - YLD_OFFSET_3 + 1; % Adjust index start (+1)
    yld_refine_end = yld_refine_start + CHJ_SIGNAL_LENGTH - 1;

    % Add boundary checks for YLD segment extraction
    if yld_refine_start < 1 || yld_refine_end > len_yld1
         fprintf('Event %d: Warning! YLD refinement window [%d, %d] out of bounds. Skipping event.\n', i, yld_refine_start, yld_refine_end);
         continue;
    end
    processed_yld_signal = filtered_yld_signal1(yld_refine_start : yld_refine_end);
    % Apply detrending and windowing (assuming windowsignal function exists)
    if exist('windowsignal', 'file')
       processed_yld_signal = real(windowsignal(detrend(processed_yld_signal)));
    else
       processed_yld_signal = real(detrend(processed_yld_signal)); % Fallback if no windowing function
       warning('Function "windowsignal" not found. Applying detrend only.');
    end

    % Define sliding window steps for refinement
    subsignal_step = floor(CHJ_SIGNAL_LENGTH / 4); % Step size (e.g., 1/4 window) - ADJUST AS NEEDED
    extracted_len = length(chj_match_signal1);
    subsignal_starts_relative = 1 : subsignal_step : (extracted_len - CHJ_SIGNAL_LENGTH + 1); % Relative indices within extracted segment

    for subi = 1:numel(subsignal_starts_relative)
        current_sub_start_relative = floor(subsignal_starts_relative(subi)); % Relative start index
        current_sub_end_relative = current_sub_start_relative + CHJ_SIGNAL_LENGTH - 1;

        % Get CHJ segments for this sub-window
        processed_chj_signal1 = chj_match_signal1(current_sub_start_relative : current_sub_end_relative);
        processed_chj_signal2 = chj_match_signal2(current_sub_start_relative : current_sub_end_relative);
        processed_chj_signal3 = chj_match_signal3(current_sub_start_relative : current_sub_end_relative);

        % Preprocess CHJ segments
        if exist('windowsignal', 'file')
          processed_chj_signal1 = real(windowsignal(detrend(processed_chj_signal1)));
          processed_chj_signal2 = real(windowsignal(detrend(processed_chj_signal2)));
          processed_chj_signal3 = real(windowsignal(detrend(processed_chj_signal3)));
        else
           processed_chj_signal1 = real(detrend(processed_chj_signal1)); % Fallback
           processed_chj_signal2 = real(detrend(processed_chj_signal2));
           processed_chj_signal3 = real(detrend(processed_chj_signal3));
        end

        % Refined Cross-Correlation
        [r_gcc, lags_gcc] = xcorr(processed_chj_signal1, processed_yld_signal, 'normalized');
        [R_gcc, max_lag_idx] = max(r_gcc);

        % Get lag in samples (assuming cal_tau provides this, or use lags_gcc)
        if exist('cal_tau', 'file')
            t_gcc_samples = cal_tau(r_gcc, lags_gcc'); % Use your function
        else
            t_gcc_samples = lags_gcc(max_lag_idx); % Lag at max correlation (samples)
            warning('Function "cal_tau" not found. Using simple lag at max R.');
        end

        % Check refined correlation threshold
        if R_gcc < REFINE_CORR_THR
            continue
        end

        % Get CHJ Direction of Arrival for this sub-window
        % Assuming get_2d_result_single_window exists and returns angles in degrees
        try
            % Pass absolute start index of the coarse window for potential timing reference inside function
            [~, chj_azimuth, chj_elevation, ~, ~] = get_2d_result_single_window(coarse_chj_loc, processed_chj_signal1, processed_chj_signal2, processed_chj_signal3);
        catch ME_DOA
            fprintf('Event %d: Error in get_2d_result_single_window for sub-window %d. Skipping sub-window.\n', i, subi);
            fprintf('Error message: %s\n', ME_DOA.message);
            continue;
        end

        % Check if DOA calculation returned valid angles (e.g., not NaN or placeholder)
        if isnan(chj_azimuth) || isnan(chj_elevation) % Adjust condition based on function's failure output
            fprintf('Event %d: Invalid DOA results for sub-window %d. Skipping.\n', i, subi);
            continue
        end

        % --- D. Geometric 3D Localization ---
        % Convert angles to radians and then to Cartesian unit vectors
        % Verify angle conventions: sph2cart(AZ, EL, R)
        % AZ=azimuth angle, EL=elevation angle. Check if 0 degrees elevation is horizon or zenith.
        % Assuming 0 elevation = horizon, 0 azimuth = X-axis? Adjust if needed.
        % Your original code used 90-elevation, implying elevation from zenith? Using that here.
        [R1_x, R1_y, R1_z] = sph2cart(deg2rad(yld_azimuth(i)), deg2rad(90-yld_elevation(i)), 1); % YLD direction
        [R2_x, R2_y, R2_z] = sph2cart(deg2rad(chj_azimuth), deg2rad(90-chj_elevation), 1); % CHJ direction

        A1 = [R1_x; R1_y; R1_z]; % YLD Direction vector (Column)
        A2 = [R2_x; R2_y; R2_z]; % CHJ Direction vector (Column)

        % Matrix for solving system M * [R1; -R2; R3] = p
        C = cross(A1, A2);
        normC = norm(C);
        if normC < 1e-9 % Check for nearly parallel vectors (increased tolerance slightly)
             % fprintf('Event %d: Warning! Direction vectors nearly parallel for sub-window %d. Skipping.\n', i, subi);
             continue;
        end
        c_unit = C / normC;
        M = [A1, -A2, c_unit];

        % Check conditioning before solving
        if rcond(M) < 1e-12 % Check reciprocal condition number (increased tolerance slightly)
             % fprintf('Event %d: Warning! Matrix M ill-conditioned for sub-window %d (rcond=%.2e). Skipping.\n', i, subi, rcond(M));
             continue;
        end

        % Solve for scalars using backslash (more stable)
        try
            sol = M \ P_BASELINE(:); % Ensure baseline is column vector
            R1_value = sol(1);      % Distance along YLD ray
            R2_value = -sol(2);     % Distance along CHJ ray (note sign change from M)
        catch ME_SOLVE
            fprintf('Event %d: Error solving linear system M\\p for sub-window %d. Skipping.\n', i, subi);
            fprintf('Error message: %s\n', ME_SOLVE.message);
            continue;
        end

        % Physical check: Distances from stations should be positive
        if R1_value < 0 || R2_value < 0
             % fprintf('Event %d: Warning! Negative distance calculated for sub-window %d (R1=%.1f, R2=%.1f). Skipping.\n', i, subi, R1_value, R2_value);
             continue;
        end

        % Calculate source location S as midpoint of the segment connecting the closest points on rays
        closest_point_on_ray1 = YLD_SIT(:) + R1_value * A1;
        closest_point_on_ray2 = CHJ_SIT(:) + R2_value * A2;
        sub_S = (closest_point_on_ray1 + closest_point_on_ray2) / 2;

        % --- E. Time Consistency Validation ---
        t_yld_geom = norm(sub_S - YLD_SIT(:)) / C_LIGHT; % Propagation time from S to YLD (seconds)
        t_chj_geom = norm(sub_S - CHJ_SIT(:)) / C_LIGHT; % Propagation time from S to CHJ (seconds)
        dlta_t_geom_s = t_yld_geom - t_chj_geom; % Theoretical TDOA based on geometry (seconds) - Note: YLD - CHJ

        % Measured TDOA from correlation lag (convert samples to ns)
        % Assumes t_gcc_samples is lag of CHJ relative to YLD (CHJ time - YLD time)
        dlta_T_measured_ns = t_gcc_samples * TGCC_TIME_MULTIPLIER;

        % Convert geometric TDOA to ns for comparison
        dlta_t_geom_ns = dlta_t_geom_s * 1e9;

        % Calculate residual (absolute difference) in nanoseconds
        dlta_residual_ns = abs(dlta_t_geom_ns - dlta_T_measured_ns);

        dltas_ns = [dltas_ns; dlta_residual_ns]; % Store residual for analysis

        % Check if residual is within tolerance
        if dlta_residual_ns <= TIME_VALIDATION_TOLERANCE_W
            % Store valid results for this sub-window
            sub_S_results = [sub_S_results; sub_S']; % Store location (as row)
            sub_R_gccs = [sub_R_gccs; R_gcc];       % Store correlation
            sub_refined_indices = [sub_refined_indices; refine_search_start + current_sub_start_relative - 1]; % Store absolute start index of refined window
            sub_t_gcc_samples = [sub_t_gcc_samples; t_gcc_samples]; % Store lag
        end
    end % End of sliding window refinement loop (subi)

    % --- F. Select Best Result from Refinement ---
    if ~isempty(sub_R_gccs)
        [max_R_gcc_refined, best_idx] = max(sub_R_gccs);
        best_S = sub_S_results(best_idx, :);
        best_refined_chj_start_index = sub_refined_indices(best_idx);
        best_time_residual = abs( (norm(best_S - YLD_SIT) - norm(best_S - CHJ_SIT))/C_LIGHT*1e9 - sub_t_gcc_samples(best_idx)*TGCC_TIME_MULTIPLIER );

        S_results = [S_results; best_S]; % Store the best location for this event i

        % Store detailed results for the best sub-window
        match_results(end+1).yld_start_loc = yld_start_loc(i);
        match_results(end).chj_loc_refined_start = best_refined_chj_start_index; % Absolute index of refined window start
        match_results(end).r_gccs_refined = max_R_gcc_refined;
        match_results(end).S_calc = best_S;
        match_results(end).time_residual_ns = best_time_residual;

        fprintf('Event %d: Successfully processed. Best refined R=%.3f. Stored location: [%.1f, %.1f, %.1f]. Time residual: %.1f ns\n', ...
                i, max_R_gcc_refined, best_S(1), best_S(2), best_S(3), best_time_residual);
    else
        fprintf('Event %d: No valid 3D location found meeting criteria after refinement.\n', i);
    end

end % End of main event loop (i)

close(h); % Close waitbar
fprintf('Processing finished. %d valid events found and localized.\n', size(S_results, 1));

%% Step 3: Post-processing & Plotting ====================================
fprintf('Step 3: Filtering results and plotting...\n');

if isempty(S_results)
    fprintf('No valid S_results to plot.\n');
    return;
end

% Apply spatial filtering to final results (optional)
x_range = [-50000, 50000]; % Example reasonable range for X (meters)
y_range = [-50000, 0];     % Example reasonable range for Y (meters)
z_range = [0, 50000];      % Example reasonable range for Z (meters, altitude > 0)

filter_mask = S_results(:,1) >= x_range(1) & S_results(:,1) <= x_range(2) & ...
              S_results(:,2) >= y_range(1) & S_results(:,2) <= y_range(2) & ...
              S_results(:,3) >= z_range(1) & S_results(:,3) <= z_range(2);

filtered_S = S_results(filter_mask, :);
fprintf('Applied spatial filter: %d out of %d points remain.\n', size(filtered_S, 1), size(S_results, 1));

if isempty(filtered_S)
    fprintf('No points remaining after spatial filter.\n');
    return;
end

% --- Plotting ---
% 1. 3D Scatter Plot
figure('Name', '3D Localization Results');
scatter3(filtered_S(:,1), filtered_S(:,2), filtered_S(:,3), 10, filtered_S(:,3), 'filled'); % Color by altitude (Z)
hold on;
scatter3(YLD_SIT(1), YLD_SIT(2), YLD_SIT(3), 100, 'r^', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'YLD Station');
scatter3(CHJ_SIT(1), CHJ_SIT(2), CHJ_SIT(3), 100, 'bs', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'CHJ Station');
hold off;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title(sprintf('3D Lightning Locations (%d points)', size(filtered_S, 1)));
colorbar; % Shows altitude scale
axis equal; grid on; view(3); % Set 3D view
legend('show', 'Location', 'best');

% 2. Azimuth-Elevation Plot (relative to YLD station)
figure('Name', 'Azimuth-Elevation Plot (from YLD)');
% Calculate vectors from YLD station to source points
vectors_from_yld = filtered_S - YLD_SIT;
[azimuth, elevation, ~] = cart2sph(vectors_from_yld(:,1), vectors_from_yld(:,2), vectors_from_yld(:,3));
azimuth = mod(azimuth, 2*pi);
% 然后，将角度转换为以正北方向为参考的顺时针角度
% 正北方向对应的角度是pi/2
azimuth_north = pi/2 - azimuth;
% 由于azimuth_north的范围可能超出[-pi, pi]，需要进行调整
azimuth_north = mod(azimuth_north + pi, 2*pi) - pi;
% 将弧度转换为角度
azimuth = azimuth_north * 180 / pi;
elevation = elevation * 180 / pi;
azimuth = mod(azimuth, 360);
% 调整仰角范围为 0-90 度
elevation = abs(elevation); % 确保仰角为非负值
elevation(elevation > 90) = 90; % 限制最大值为 90 度

scatter(azimuth, elevation, 10, filtered_S(:,3), 'filled'); % Color by altitude
xlabel('Azimuth (degrees) from YLD');
ylabel('Elevation (degrees) from YLD');
title('Source Azimuth/Elevation relative to YLD Station');
grid on;
colorbar; % Shows altitude scale
xlim([0 360]); % Standard azimuth range
ylim([-90 90]); % Standard elevation range

% 3. Top-Down View (X-Y Plane)
figure('Name', 'Top-Down View (X-Y Plane)');
scatter(filtered_S(:,1), filtered_S(:,2), 10, filtered_S(:,3), 'filled'); % Color by altitude
hold on;
plot(YLD_SIT(1), YLD_SIT(2), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'YLD Station');
plot(CHJ_SIT(1), CHJ_SIT(2), 'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'CHJ Station');
hold off;
xlabel('X (m)');
ylabel('Y (m)');
title('Top-Down View (X-Y)');
colorbar; % Shows altitude scale
axis equal; grid on;
legend('show', 'Location', 'best');

fprintf('Plotting complete.\n');

% Optional: Save results
% save('localization_results.mat', 'S_results', 'match_results', 'dltas_ns');