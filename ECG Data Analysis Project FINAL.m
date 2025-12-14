%% WEEK 15 PROJECT - ECG ANALYSIS
% Assumes the following files are in the current folder:
%   Resting.mat, Exercise.mat, BoxBreathing.mat
% Each file must contain a variable b1 with columns:
%   1: time (s)
%   2: ECG (mV)
%   3: instantaneous heart rate (not used directly)

clear; clc; close all;

%% STEP 1: LOAD DATA
% Load Resting
load('Resting.mat');      % loads b1
t_rest   = b1(:,1);
ecg_rest = b1(:,2);
hr_rest_raw = b1(:,3);    %#ok<NASGU> % not used, but kept for reference
clear b1;

% Load Exercise
load('Exercise.mat');
t_ex   = b1(:,1);
ecg_ex = b1(:,2);
hr_ex_raw = b1(:,3);      %#ok<NASGU>
clear b1;

% Load Box Breathing
load('BoxBreathing.mat');
t_box   = b1(:,1);
ecg_box = b1(:,2);
hr_box_raw = b1(:,3);     %#ok<NASGU>
clear b1;

% Sampling frequencies (approx.)
fs_rest = 1/mean(diff(t_rest));
fs_ex   = 1/mean(diff(t_ex));
fs_box  = 1/mean(diff(t_box));

%% STEP 2: PLOT 10 SECONDS OF CLEAN ECG FOR EACH CONDITION
% Just use the first 10 seconds

% Resting - 10s
t0 = t_rest(1);
idx10_rest = t_rest >= t0 & t_rest < t0 + 10;
figure;
plot(t_rest(idx10_rest), ecg_rest(idx10_rest), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('ECG (mV)');
title('Resting ECG - First 10 Seconds');
grid on;

% Exercise - 10s
t0 = t_ex(1);
idx10_ex = t_ex >= t0 & t_ex < t0 + 10;
figure;
plot(t_ex(idx10_ex), ecg_ex(idx10_ex), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('ECG (mV)');
title('Exercise ECG - First 10 Seconds');
grid on;

% Box Breathing - 10s
t0 = t_box(1);
idx10_box = t_box >= t0 & t_box < t0 + 10;
figure;
plot(t_box(idx10_box), ecg_box(idx10_box), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('ECG (mV)');
title('Box Breathing ECG - First 10 Seconds');
grid on;

%% STEP 2B: 3x1 SUBPLOT WITH ~60 SECONDS OF EACH
dur60 = 60;  % seconds

idx60_rest = t_rest >= t_rest(1) & t_rest < t_rest(1) + dur60;
idx60_ex   = t_ex   >= t_ex(1)   & t_ex   < t_ex(1)   + dur60;
idx60_box  = t_box  >= t_box(1)  & t_box  < t_box(1)  + dur60;

figure;
subplot(3,1,1);
plot(t_rest(idx60_rest), ecg_rest(idx60_rest), 'LineWidth', 1);
ylabel('ECG (mV)');
title('Resting (~60 s)');
grid on;

subplot(3,1,2);
plot(t_ex(idx60_ex), ecg_ex(idx60_ex), 'LineWidth', 1);
ylabel('ECG (mV)');
title('Exercise (~60 s)');
grid on;

subplot(3,1,3);
plot(t_box(idx60_box), ecg_box(idx60_box), 'LineWidth', 1);
xlabel('Time (s)');
ylabel('ECG (mV)');
title('Box Breathing (~60 s)');
grid on;

%% STEP 3 & 4: PEAK DETECTION OVERVIEW USING findpeaks (DEFAULT SETTINGS)
% Use small Resting window (5 seconds) for demonstration

demo_start = t_rest(1) + 20;      % start time in seconds (adjust if needed)
demo_end   = demo_start + 5;      % 5 second window
demo_idx   = t_rest >= demo_start & t_rest <= demo_end;

t_demo   = t_rest(demo_idx);
ecg_demo = ecg_rest(demo_idx);

figure;
subplot(2,1,1);
plot(t_demo, ecg_demo, 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('ECG (mV)');
title('Resting ECG - Demo Window');
grid on;

% Default findpeaks (no special parameters)
[pks_demo, locs_demo] = findpeaks(ecg_demo, t_demo);

subplot(2,1,2);
plot(t_demo, ecg_demo, 'LineWidth', 1.2); hold on;
plot(locs_demo, pks_demo, 'ro', 'MarkerSize', 6, 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('ECG (mV)');
title('findpeaks() with Basic Settings - Detects P, QRS, T etc.');
grid on;

%% STEP 5: DETECT ONLY R PEAKS (CHOOSE APPROPRIATE PARAMETERS)
% Height-based threshold, scaled by the max amplitude

minHeight_rest = 0.5 * max(ecg_rest);   % adjust if needed
minHeight_ex   = 0.5 * max(ecg_ex);
minHeight_box  = 0.5 * max(ecg_box);

% Minimum distance between beats (in seconds) -> convert to samples
minRR_sec  = 0.3;   % 0.3 s ~ 200 bpm, should work for both rest & exercise
minRR_rest = round(minRR_sec * fs_rest);
minRR_ex   = round(minRR_sec * fs_ex);
minRR_box  = round(minRR_sec * fs_box);

% Use sample-based findpeaks to get indices of R peaks
[Rpks_rest, Rlocs_rest] = findpeaks(ecg_rest, ...
    'MinPeakHeight', minHeight_rest, ...
    'MinPeakDistance', minRR_rest);

[Rpks_ex, Rlocs_ex] = findpeaks(ecg_ex, ...
    'MinPeakHeight', minHeight_ex, ...
    'MinPeakDistance', minRR_ex);

[Rpks_box, Rlocs_box] = findpeaks(ecg_box, ...
    'MinPeakHeight', minHeight_box, ...
    'MinPeakDistance', minRR_box);

% Show R peaks on a small Resting window
figure;
plot(t_rest, ecg_rest, 'LineWidth', 1); hold on;
plot(t_rest(Rlocs_rest), Rpks_rest, 'ro', 'MarkerSize', 4);
xlim([demo_start demo_end]);   % zoom near earlier demo window
xlabel('Time (s)');
ylabel('ECG (mV)');
title('Resting ECG with Only R Peaks Detected');
grid on;

%% STEP 6: ALIGN HEARTBEATS AROUND R PEAKS (REST, EXERCISE, BOX)
% Fixed window around each R peak (in samples)
win_before = round(0.25 * fs_rest);   % 250 ms before
win_after  = round(0.45 * fs_rest);   % 450 ms after

% Helper function (anonymous) to collect aligned beats
collect_beats = @(ecg, Rlocs, fs, wb, wa) ...
    local_collect_beats(ecg, Rlocs, fs, wb, wa);

[beats_rest, t_rel_rest] = collect_beats(ecg_rest, Rlocs_rest, fs_rest, win_before, win_after);
[beats_ex,   t_rel_ex]   = collect_beats(ecg_ex,   Rlocs_ex,   fs_ex,   win_before, win_after);
[beats_box,  t_rel_box]  = collect_beats(ecg_box,  Rlocs_box,  fs_box,  win_before, win_after);

% Plot aligned heartbeats for each condition
figure;
subplot(3,1,1);
plot(t_rel_rest, beats_rest, 'LineWidth', 0.5);
xlabel('Time relative to R peak (s)');
ylabel('ECG (mV)');
title('Aligned Heartbeats - Resting');
grid on;

subplot(3,1,2);
plot(t_rel_ex, beats_ex, 'LineWidth', 0.5);
xlabel('Time relative to R peak (s)');
ylabel('ECG (mV)');
title('Aligned Heartbeats - Exercise');
grid on;

subplot(3,1,3);
plot(t_rel_box, beats_box, 'LineWidth', 0.5);
xlabel('Time relative to R peak (s)');
ylabel('ECG (mV)');
title('Aligned Heartbeats - Box Breathing');
grid on;

%% STEP 7: HEART RATE VARIABILITY (HRV) FROM R-R INTERVALS
% Use a specific 10-second window (120 to 130)

t_start = 120;
t_end   = 130;

% R-peak times
Rtimes_rest = t_rest(Rlocs_rest);
Rtimes_ex   = t_ex(Rlocs_ex);
Rtimes_box  = t_box(Rlocs_box);

% Logical index for peaks in that window
winR_rest = Rtimes_rest >= t_start & Rtimes_rest <= t_end;
winR_ex   = Rtimes_ex   >= t_start & Rtimes_ex   <= t_end;
winR_box  = Rtimes_box  >= t_start & Rtimes_box  <= t_end;

Rwin_rest = Rtimes_rest(winR_rest);
Rwin_ex   = Rtimes_ex(winR_ex);
Rwin_box  = Rtimes_box(winR_box);

% Interbeat intervals (IBI) - manual subtraction
IBI_rest_manual = Rwin_rest(2:end) - Rwin_rest(1:end-1);
IBI_ex_manual   = Rwin_ex(2:end)   - Rwin_ex(1:end-1);
IBI_box_manual  = Rwin_box(2:end)  - Rwin_box(1:end-1);

% Using diff() - should match manual above
IBI_rest_diff = diff(Rwin_rest);
IBI_ex_diff   = diff(Rwin_ex);
IBI_box_diff  = diff(Rwin_box);

% Convert intervals (seconds) to beats per minute
HR_rest = 60 ./ IBI_rest_diff;
HR_ex   = 60 ./ IBI_ex_diff;
HR_box  = 60 ./ IBI_box_diff;

% Histograms for heart rate in each condition
figure;
subplot(3,1,1);
histogram(HR_rest);
xlabel('Heart Rate (bpm)');
ylabel('Count');
title('Heart Rate Histogram - Resting');

subplot(3,1,2);
histogram(HR_ex);
xlabel('Heart Rate (bpm)');
ylabel('Count');
title('Heart Rate Histogram - Exercise');

subplot(3,1,3);
histogram(HR_box);
xlabel('Heart Rate (bpm)');
ylabel('Count');
title('Heart Rate Histogram - Box Breathing');

% Boxchart comparing Resting vs Exercise
group = [repmat("Rest", size(HR_rest)); repmat("Exercise", size(HR_ex))];
HR_all = [HR_rest(:); HR_ex(:)];

figure;
boxchart(group, HR_all);
ylabel('Heart Rate (bpm)');
title('Heart Rate Comparison: Resting vs Exercise');

%% STEP 8: SLOPES AROUND R AND T PULSES (RESTING ONLY)
% Very simple approximation:
%   R peaks: we already have exact indices (Rlocs_rest)
%   T peaks: approximate as a fixed time (~200 ms) after each R peak

offset_T_sec = 0.2;                        % 200 ms after R peak
offset_T_samp = round(offset_T_sec * fs_rest);

% Use 10 samples before and after each point for slope calculation
N = 10;   % samples before and after

numR = length(Rlocs_rest);
slopes_R = [];
slopes_T = [];

for k = 1:numR
    idxR = Rlocs_rest(k);
    idxT = idxR + offset_T_samp;   % approximate T location

    % R slope (only if window is inside data)
    if idxR - N >= 1 && idxR + N <= length(ecg_rest)
        dt_R = (2*N) / fs_rest;  % total window time (s)
        dECG_R = ecg_rest(idxR + N) - ecg_rest(idxR - N);
        slopes_R(end+1) = abs(dECG_R / dt_R); %#ok<SAGROW>
    end

    % T slope (only if window is inside data)
    if idxT - N >= 1 && idxT + N <= length(ecg_rest)
        dt_T = (2*N) / fs_rest;
        dECG_T = ecg_rest(idxT + N) - ecg_rest(idxT - N);
        slopes_T(end+1) = abs(dECG_T / dt_T); %#ok<SAGROW>
    end
end

figure;
histogram(slopes_R); hold on;
histogram(slopes_T);
xlabel('Absolute Slope (mV/s)');
ylabel('Count');
title('Slope Around R vs T Pulses (Resting)');
legend('R slopes','T slopes');
grid on;

%% STEP 9: STATISTICS WITH ttest2
% Test 1: Rest vs Exercise HR (same 120-130s window)

[h1, p1] = ttest2(HR_rest, HR_ex);

fprintf('\nT-TEST 1: Rest vs Exercise (HR in %.0f-%.0f s window)\n', t_start, t_end);
fprintf('  h = %d (1 = reject null, 0 = fail to reject)\n', h1);
fprintf('  p-value = %.4f\n', p1);

if h1 == 1
    fprintf('  => Rest and Exercise heart rates ARE statistically different.\n');
else
    fprintf('  => Rest and Exercise heart rates are NOT statistically different.\n');
end

% Test 2: Two different time windows within RESTING (should be similar)
% Window A: same as before (t_start to t_end)
% Window B: 10 seconds later (t_start+20 to t_end+20)

t_start_B = t_start + 20;
t_end_B   = t_end   + 20;

winR_rest_B = Rtimes_rest >= t_start_B & Rtimes_rest <= t_end_B;
Rwin_rest_B = Rtimes_rest(winR_rest_B);
IBI_rest_B  = diff(Rwin_rest_B);
HR_rest_B   = 60 ./ IBI_rest_B;

% Make sure we actually have data in both windows
if ~isempty(HR_rest) && ~isempty(HR_rest_B)
    [h2, p2] = ttest2(HR_rest, HR_rest_B);

    fprintf('\nT-TEST 2: Rest Window A vs Rest Window B\n');
    fprintf('  Window A: %.0f-%.0f s,  Window B: %.0f-%.0f s\n', ...
        t_start, t_end, t_start_B, t_end_B);
    fprintf('  h = %d (1 = reject null, 0 = fail to reject)\n', h2);
    fprintf('  p-value = %.4f\n', p2);

    if h2 == 1
        fprintf('  => The two resting windows ARE statistically different.\n');
    else
        fprintf('  => The two resting windows are NOT statistically different.\n');
    end
else
    fprintf('\nT-TEST 2: Not enough data in one or both resting windows.\n');
end

%% LOCAL FUNCTION (AT END OF SCRIPT)
% Simple helper collects aligned beats around each R peak
function [beats, t_rel] = local_collect_beats(ecg, Rlocs, fs, win_before, win_after)
    nSamples = length(ecg);
    beatLen  = win_before + win_after + 1;
    beats    = [];

    for i = 1:length(Rlocs)
        idx = Rlocs(i);
        if idx - win_before >= 1 && idx + win_after <= nSamples
            seg = ecg(idx - win_before : idx + win_after);
            beats = [beats; seg']; %#ok<AGROW>
        end
    end

    t_rel = (-win_before:win_after) / fs;
    beats = beats';  % each column is one beat
end