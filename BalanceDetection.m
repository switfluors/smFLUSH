%==========================================================================
% BM_G1_balanced_v2_abc.m
%
% Single-molecule spectral-centroid analysis (Balanced Method).
%
% Pipeline:
%   1. Load a multi-frame ND2 stack (Bio-Formats) and build a wavelength
%      calibration from known reference lines.
%   2. Locate the spatial slit and estimate / subtract a global background.
%   3. Detect single-molecule events (spatial peaks per frame) and extract
%      their dispersed spectra onto a common wavelength grid.
%   4. Filter events by photon budget and refine the spectral centroid with
%      the iterative Balanced Method (fixed-width window centered on the
%      previous centroid; window half-width set once by the 1/N method).
%   5. Produce a 3-panel summary figure (a: spectra heatmap, b: mean
%      spectrum, c: centroid histogram) and log results to .mat, .csv and a
%      multi-sheet Excel workbook.
%
% Requirements: Bio-Formats toolbox (bfmatlab) on the MATLAB path.
%==========================================================================

%% Spectral calibration and loading
clear;
clc;

% Wavelength calibration: known emission lines (xi, nm) at measured detector
% rows (yi, px). A linear fit maps wavelength -> row coordinate.
yi = [176, 190, 206, 225];          % detector row (px) of each reference line
xi = [436.6, 487.7, 546.5, 611.6];  % corresponding wavelength (nm)
x0 = 53;                            % horizontal offset (px) of the 0th order

fy  = 400:3:900;   % coarse wavelength grid (extended floor 500->400 nm to capture BODIPY-517 emission onset)
fy2 = 400:1:900;   % fine (1 nm) wavelength grid used for all spectral analysis
coeffs = polyfit(xi, yi, 1);        % linear wavelength -> row calibration
fx  = polyval(coeffs, fy);          % row coordinate for the coarse grid
fx2 = polyval(coeffs, fy2);         % row coordinate for the fine grid

addpath('D:\bfmatlab')

% Input ND2 stack.
nd2FilePath = 'D:\Miami\MIN-BOP_B6\G1\G1_006.nd2';

% Derive the output base name from the actual input file, so every output
% (.mat, .pdf, .csv, Excel rows) is tagged after the file being analyzed.
[~, fileName, ext] = fileparts(nd2FilePath);
fileNameWithEscapes = strrep(fileName, '_', '__');   % escape '_' for on-plot text (avoids TeX subscript)
fullFileName = [fileNameWithEscapes, ext];

% Output folder, shared by the .mat, .pdf, .csv and Excel log below.
outputFilePath = 'D:\Miami\MIN-BOP_B6\G1\G1.xlsx';   % full Excel path; folder is derived from it
[outputFolder, ~, ~] = fileparts(outputFilePath);
if ~isempty(outputFolder) && ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Load the ND2 stack via Bio-Formats.
data = bfopen(nd2FilePath);
numFrames = size(data{1}, 1);                       % total number of frames
[frameHeight, frameWidth] = size(data{1}{1});       % dimensions of a single frame

% Assemble all frames into a 3-D array (height x width x frame).
A = zeros(frameHeight, frameWidth, numFrames);
for i = 1:numFrames
    A(:, :, i) = data{1}{i};
end
A = rot90(A, 4);   % orientation correction (4x90 deg = identity; kept for parity with acquisition convention)
A9 = A;

%% Slit localization and preview
% Average image across all frames.
Ave1 = mean(A9, 3);

% Locate the slit: in each row, the brightest pixel within the first 100
% columns marks the slit; the mode over all rows is the slit column.
slit_position = [];
for n = 1:size(Ave1, 1)
    [m, I] = max(Ave1(n, 1:100));
    slit_position = [slit_position, I];
end
slit_position = slit_position';
slit = mode(slit_position);

% Preview: average image (left) and a short frame animation (right).
figure;
subplot(1, 2, 1);
imagesc(Ave1);
colormap(gray);

subplot(1, 2, 2);
for n = 1:200
    imagesc(A(:, :, n));
    title(['Frame ', num2str(n)]);
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    colormap(gray);
    set(gca, 'color', 'none', 'Position', [0.55 0.1 0.4 0.8]);
    drawnow;
    % pause(0.02);
end

%% Global background estimation
% For each row, sample columns 150-350 over all frames; frames whose row-mean
% is below threshold (504) are treated as signal-free and averaged to form
% the per-row background profile.
back = [];
for n = 1:size(A, 1)
    w = squeeze(A(n, 150:350, 1:end));
    w = mean(w, 1);
    a = w < 504;                       % signal-free frame mask for this row
    back(n, :) = mean(A(n, :, a), 3);
end
clf
imagesc(back)
colormap(gray)

%% Background subtraction
back = repmat(back, [1, 1, size(A, 3)]);
A = A - back;

%% Peak (event) detection
% In each frame, sum a narrow band around the slit and find spatial peaks;
% for each peak, locate its dispersion maximum to record (row, col, frame).
TR = [];
for n = 1:size(A, 3)
    b = A(:, slit-3:slit+3, n);
    s = sum(b, 2);
    t = 800;                           % minimum peak height (photon threshold)
    warning off
    [~, y] = findpeaks(s, 'MinPeakHeight', t, 'MinPeakDistance', 3);
    if size(y, 2) >= 1
        for m = 1:size(y, 1)
            if y < 196 & y > 3         % keep events away from the frame edges
                b2 = A(y(m)-2:y(m)+2, 1:slit+10, n);
                s_x = sum(b2, 1);
                [~, x] = max(s_x);
                IF = [y(m), x, n];
                TR = [TR; IF];
            end
        end
    end
end

%% Spectral extraction onto the common wavelength grid
s = TR(:, [2, 1, 3]);                  % reorder to (col, row, frame)
y0 = 0;
bb3 = [];
for n = 1:size(s, 1)
    % Extract a 3-row band along the dispersion axis for this event.
    bb = A((-1:1)+s(n,2)+y0, (round(fx(1)):round(fx(end)))-x0+s(n,1), s(n,3));
    bb2 = sum(bb, 1);                  % collapse to a 1-D spectrum
    bb3 = [bb3, sum(bb2)];             % integrated intensity per event
    bs(n, :) = interp1([round(fx(1)):round(fx(end))], bb2, fx2);  % resample onto fine grid
end

%% Balanced Method: iterative spectral-centroid refinement
num          = 10;      % maximum number of refinement rounds
convTol      = 0.0001;  % nm; stop early once consecutive rounds agree within this tolerance
level_frac   = 1/10;    % 1/N method: window width set by the 1/N-of-peak crossing points
width_factor = 1.0;     % halfWidth = width_factor * (1/N width) / 2

bs3 = bs;
bs3(isnan(bs3)) = 0;

% Per-event photon budget (calibrated to the combined dataset).
s = bb3/100*14.4/0.95;

% Filter events by photon budget.
photonMin = 500;        % lower photon bound
photonMax = 3000;       % upper photon bound
a1 = s >= photonMin;
a2 = s <= photonMax;
a  = logical(a1 .* a2);
s   = s(a);
spe = bs3(a, :);

C        = nan(1, num); % per-round centroid estimate
numValid = nan(1, num); % per-round count of events with a usable centroid
halfWidth = nan;        % fixed window half-width, set from the 1/N method in round 1

for k = 1:num
    if k == 1
        A = 650; B = 250;          % round 1: wide window covering the whole captured range (400-900 nm)
        spe1 = spe;
    else
        A = C(k-1); B = halfWidth;  % later rounds: window centered on previous centroid, fixed half-width
        spe1 = spe - E;
    end

    % Integer column indices into spe1/fy2 for the window [A-B, A+B].
    % Indices and wavelength weights are taken from the same rounded grid so
    % they stay aligned (a fractional-index mismatch previously caused a slow
    % ~0.09 nm/round drift instead of true convergence).
    colStart = round(A - B - fy2(1) + 1);
    colEnd   = round(A + B - fy2(1) + 1);
    colStart = max(colStart, 1);
    colEnd   = min(colEnd, numel(fy2));
    wl = fy2(colStart:colEnd);          % wavelength grid for this window, aligned with spe1's columns

    % Intensity-weighted centroid per event over the window.
    c = nan(1, size(spe1, 1));
    for n = 1:size(spe1, 1)
        seg = spe1(n, colStart:colEnd);
        denom = sum(seg);
        if denom > 0
            c(n) = sum(wl.*seg) / denom;
        end
    end

    if k == 1
        % Red-edge background sample (same column as the original pipeline).
        y = mean(spe1, 1);
        E = y(numel(y) - 2);

        % --- Window half-width via the 1/N method (N = 1/level_frac) ---
        % Derived once from the round-1 averaged spectrum, then held fixed.
        y_bg = y - E;
        [peak_val, peak_idx] = max(y_bg);
        peak_x = fy2(peak_idx);         % peak wavelength (nm)
        peak_y = peak_val;              % peak intensity (a.u.)
        level = peak_val * level_frac;

        % Require nConsec consecutive below-level points before calling an
        % edge, so a single noisy background pixel cannot drag the edge out.
        nConsec = 5;
        runMask = movmin(double(y_bg < level), nConsec, 'Endpoints', 'discard');
        % runMask(i) == 1 means y_bg(i:i+nConsec-1) are ALL below level.

        % Left edge: from the peak going left, first below-level run; its
        % rightmost point (closest to the peak) is the edge (linearly interpolated).
        left_search_end = max(peak_idx - nConsec, 1);
        wl_left = fy2(1);
        for ii = left_search_end:-1:1
            if runMask(ii) == 1
                edge_idx = ii + nConsec - 1;
                if edge_idx < numel(fy2)
                    x1 = fy2(edge_idx);   y1v = y_bg(edge_idx);
                    x2 = fy2(edge_idx+1); y2v = y_bg(edge_idx+1);
                    wl_left = x1 + (level - y1v) * (x2 - x1) / (y2v - y1v);
                else
                    wl_left = fy2(edge_idx);
                end
                break;
            end
        end

        % Right edge: from the peak going right, first below-level run; its
        % leftmost point (closest to the peak) is the edge (linearly interpolated).
        right_search_start = min(peak_idx + nConsec, numel(runMask));
        wl_right = fy2(end);
        for ii = right_search_start:numel(runMask)
            if runMask(ii) == 1
                edge_idx = ii;
                if edge_idx > 1
                    x1 = fy2(edge_idx-1); y1v = y_bg(edge_idx-1);
                    x2 = fy2(edge_idx);   y2v = y_bg(edge_idx);
                    wl_right = x1 + (level - y1v) * (x2 - x1) / (y2v - y1v);
                else
                    wl_right = fy2(edge_idx);
                end
                break;
            end
        end

        width_1_100 = wl_right - wl_left;
        halfWidth   = width_factor * width_1_100 / 2;

        fprintf('\n[Balanced Method window half-width (1/%g method, from round-1 averaged spectrum)]\n', 1/level_frac);
        fprintf('  Peak = %.4f at %.2f nm\n', peak_val, fy2(peak_idx));
        fprintf('  Level (1/%g peak) = %.4f\n', 1/level_frac, level);
        fprintf('  Left edge = %.2f nm, Right edge = %.2f nm (requires %d consecutive points below level)\n', wl_left, wl_right, nConsec);
        fprintf('  Width (1/%g) = %.2f nm -> halfWidth = %.2f nm (factor=%.2f)\n\n', ...
            1/level_frac, width_1_100, halfWidth, width_factor);
    end

    C(k) = mean(c, 'omitnan');
    numValid(k) = sum(~isnan(c));

    if k > 1 && abs(C(k) - C(k-1)) < convTol
        fprintf('Centroid converged: rounds %d and %d differ by < %.3f nm; stopping early.\n', k-1, k, convTol);
        break;
    end
end

fprintf('Balanced Method stopped at round %d (max %d). Per-round centroids: %s\n', k, num, mat2str(round(C(1:k), 4)));

% Window range as a display string.
wlString = ['[', num2str(wl(1)), ':', num2str(wl(2) - wl(1)), ':', num2str(wl(end)), ']'];

%% Summary figure: a) spectra heatmap | b) mean spectrum | c) centroid histogram
labelFontSize = 11;
eventsToShow  = 200;

figure;
set(gcf, 'Units', 'inches', 'Position', [1, 1, 13, 4.5]);

% Panel layout.
hgap = 0.06; lmar = 0.07; rmar = 0.04; bmar = 0.16; tmar = 0.08;
pw   = (1 - lmar - rmar - 2*hgap) / 3;
ph   = 1 - bmar - tmar;
posA = [lmar,                bmar, pw, ph];
posB = [lmar + pw + hgap,    bmar, pw, ph];
posC = [lmar + 2*(pw+hgap),  bmar, pw, ph];

% Trim spectra to 400-800 nm for display (fy2 = 400:1:900).
colMask  = fy2 >= 400 & fy2 <= 800;
wl_plot  = fy2(colMask);
spe_plot = spe1(:, colMask);
xTickSpec = 400:100:800;

% a) Representative single-molecule spectra (heatmap).
subplot('position', posA)
nShow   = min(eventsToShow, size(spe_plot, 1));
xTickPx = xTickSpec - 400 + 1;
imagesc(spe_plot(1:nShow, :), [100 2300])
colormap(gca, hot)
set(gca, 'XTick', xTickPx, 'XTickLabel', xTickSpec)
set(gca, 'YTick', [1, nShow])
ylabel('Events');  xlabel('Wavelength (nm)')
set(gca, 'TickLength', [0.02 0.01], 'FontSize', labelFontSize);  box on
text(-0.12, 1.05, 'a', 'Units','normalized', 'FontSize',12, ...
     'FontWeight','bold', 'VerticalAlignment','top', 'Color','k')
numEvents = numel(c);

% b) Averaged single-molecule spectrum.
subplot('position', posB)
y_mean = mean(spe_plot, 1);
plot(wl_plot, y_mean, 'b-', 'LineWidth', 1.2)
xlim([400 800]);  ylim([0, max(y_mean)*1.15])
set(gca, 'XTick', xTickSpec)
yNice = ceil(max(y_mean)/10)*10;
set(gca, 'YTick', round(linspace(0, yNice, 5)))
ylabel('I (a.u.)');  xlabel('Wavelength (nm)')
set(gca, 'TickLength', [0.02 0.01], 'FontSize', labelFontSize);  box on
text(-0.12, 1.05, 'b', 'Units','normalized', 'FontSize',12, ...
     'FontWeight','bold', 'VerticalAlignment','top')

% c) Spectral-centroid histogram.
% Set the x-axis range here (multiples of 5); leave [] for automatic limits.
cXmin = [500];
cXmax = [600];
if isempty(cXmin),  cXmin = 5 * floor(min(c)/5);  end
if isempty(cXmax),  cXmax = 5 * ceil(max(c)/5);   end
cMid    = 5 * round(mean([cXmin, cXmax]) / 5);     % middle tick (multiple of 5)
cXticks = [cXmin, cMid, cXmax];                    % exactly 3 ticks

subplot('position', posC)
binWidth = 3.5;   % nm
edges = cXmin : binWidth : cXmax;
[nb, ~] = histcounts(c, edges);
pctC = nb / sum(nb) * 100;
binCenters = (edges(1:end-1) + edges(2:end)) / 2;
bar(binCenters, pctC, 1, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none')
set(gca, 'XTick', cXticks);  xlim([cXmin cXmax])
yMaxC = ceil(max(pctC)/5)*5;
set(gca, 'YTick', 0:5:yMaxC);  ylim([0 yMaxC])
xlabel('\lambda_{SC} (nm)');  ylabel('Probability %')
set(gca, 'TickLength', [0.02 0.01], 'FontSize', labelFontSize);  box on
text(-0.12, 1.05, 'c', 'Units','normalized', 'FontSize',12, ...
     'FontWeight','bold', 'VerticalAlignment','top')

%% Summary statistics and logging
averagePhotonCount = round(mean(s));
std_y = std(c);
spectralCentroid = mean(c);

matFilePath = fullfile(outputFolder, 'Balance_Detection.mat');
pdfFilePath = fullfile(outputFolder, [fileName, '.pdf']);
csvFilePath = fullfile(outputFolder, [fileName, '_summary.csv']);

exportgraphics(gcf, pdfFilePath, 'ContentType', 'vector');

% Console summary.
disp('----------------------------------------');
disp('Summary Statistics (Balanced Method)');
disp('----------------------------------------');
fprintf('File Name:              %s\n', fullFileName);
fprintf('Number of Events:       %d\n', numEvents);
fprintf('Average Photon Count:   %.2f\n', averagePhotonCount);
fprintf('Standard Deviation:     %.2f nm\n', std_y);
fprintf('Maximum Peak:           %.2f nm\n', peak_x);
fprintf('Maximum Intensity:      %.2f nm\n', peak_y);
fprintf('Spectral Centroid:      %.2f nm\n', spectralCentroid);
fprintf('Balanced Window halfWidth (1/%g method): %.2f nm\n', 1/level_frac, halfWidth);
fprintf('Wavelength Range:       %s nm\n', wlString);
disp('----------------------------------------');

save(matFilePath, 'c', 's', 'spe', 'halfWidth');

% One-row summary table for this file: saved as its own CSV and appended to
% the running multi-file Excel log below.
summaryTable = table(...
    {fullFileName}, numEvents, averagePhotonCount, std_y, peak_x, peak_y, spectralCentroid, halfWidth, {wlString}, ...
    'VariableNames', {'File Name', 'Number of Events', 'Average Photon Count', 'Standard Deviation', 'Maximum Peak', 'Maximum Intensity', 'Spectral Centroid', 'Balanced_halfWidth_nm', 'Wavelength Range'});

writetable(summaryTable, csvFilePath);
disp(['Summary statistics saved to ', csvFilePath]);

% Append the summary row to sheet 2 of the Excel workbook (create if absent).
sheetNumber = 2;
try
    if isfile(outputFilePath)
        writetable(summaryTable, outputFilePath, 'Sheet', sheetNumber, 'WriteMode', 'append', 'WriteVariableNames', true);
    else
        writetable(summaryTable, outputFilePath, 'Sheet', sheetNumber);
    end
    disp(['Summary statistics have been written to sheet ', num2str(sheetNumber), ' of ', outputFilePath]);
catch ME
    warning('NM647:ExcelWriteFailed', ['Excel summary write failed (%s). Common causes: 1) the file is currently open in Excel; ' ...
        '2) the folder lacks write permission; 3) a previous crashed run left a stuck background Excel process ' ...
        '(end EXCEL.EXE in Task Manager and retry). The plots/analysis above this point are unaffected.'], ME.message);
end

% Per-round convergence log (one row per iteration actually run), on a
% separate sheet. Supports plots of centroid vs round or |diff| vs round.
roundNum            = (1:k)';
centroidPerRound    = C(1:k)';
diffPerRound        = [NaN; abs(diff(centroidPerRound))];
validEventsPerRound = numValid(1:k)';

iterationLog = table(repmat({fullFileName}, k, 1), roundNum, centroidPerRound, diffPerRound, validEventsPerRound, ...
    'VariableNames', {'File Name', 'Round', 'Centroid_nm', 'AbsDiff_from_prev_nm', 'Valid_Events'});

iterationSheetNumber = 3;   % separate from the one-row-per-file summary on sheet 2
try
    if isfile(outputFilePath)
        writetable(iterationLog, outputFilePath, 'Sheet', iterationSheetNumber, 'WriteMode', 'append', 'WriteVariableNames', true);
    else
        writetable(iterationLog, outputFilePath, 'Sheet', iterationSheetNumber);
    end
    disp(['Per-round convergence log written to sheet ', num2str(iterationSheetNumber), ' of ', outputFilePath]);
catch ME
    warning('NM647:IterationLogWriteFailed', 'Could not write the per-round convergence log (%s). The summary stats above are unaffected.', ME.message);
end
