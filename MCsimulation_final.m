%==========================================================================
% Monte-Carlo Simulation.m
%
% PART 1  Homogeneous emitters        -> Std_3D (Width x Psig x Center)
% PART 2  Spectral heterogeneity      -> Std_4D (Sigma x Width x Psig x Center)
%==========================================================================
warning('off', 'all');
clear; clc;
rng('shuffle');
%==========================================================================
% Experimental parameters and preloaded data
%==========================================================================
lambda = 670; NA = 1.49; px = 160/1.5; res = 2.3438;
numSpectra = 1000;
eta = 1; Nf = 1;
n1 = 16; c1 = floor(n1/2);

Psig_File_Vals     = 100:500:5000;
Center_Wavelengths = 550:5:800;
Width_Vals         = 50:50:300;
Sigma_Shift_Vals   = [0, 0.5, 1.5, 2.5];   % nm, peak-shift std (for Part 2)

%==========================================================================
% Set false to run Part 1 (homogeneous) ONLY and skip Part 2 (adding heterogeneity)
ENABLE_HETERO = true;
%==========================================================================

% Gaussian PSF
w   = 0.61*lambda/NA;  w = w/2.355;
psf = fspecial('gaussian', 16, w/px);

% Reference emission spectrum (AF647)
load 'Alexa_Fluor_647.csv';
spt_raw = Alexa_Fluor_647(:,3);
wl      = Alexa_Fluor_647(:,1);

% Wavelength axis and per-pixel background
wl_fixed  = 400:res:900;
n3        = length(wl_fixed);
n3_old    = 128;
Pbg       = 10000 * (n3 / n3_old);
Ibg       = Pbg * eta;

GT_Centroid_Fixed = sum(wl(:) .* spt_raw(:)) / sum(spt_raw(:));

font_title = 14; font_label = 12; font_tick = 10; font_cb = 10;

num_P = length(Psig_File_Vals);
num_W = length(Width_Vals);
num_C = length(Center_Wavelengths);
num_S = length(Sigma_Shift_Vals);

%==========================================================================
% PART 1 : HOMOGENEOUS EMITTERS
%==========================================================================
spt = repmat(spt_raw', numSpectra, 1);     % identical emitters
wl  = wl';                                 % row vector for conv/interp

Std_3D      = zeros(num_W, num_P, num_C);
Accuracy_3D = zeros(num_W, num_P, num_C);

disp('Part 1: one realization per Psig, sweep all windows...');
total_runs = num_P;  run_count = 0;
f1 = waitbar(0, 'Part 1...');

for i = 1:num_P
    Isig = Psig_File_Vals(i) * eta;
    sptimg4         = zeros(16, n3, numSpectra);
    First_four_rows = zeros(4,  n3, numSpectra);
    for n = 1:numSpectra
        sptimg = conv2(psf, spt(n,:));
        sptimg(:, end-(c1-2):end) = []; sptimg(:, 1:c1) = [];
        sptimg2 = interp1(wl, sptimg', wl_fixed)';
        sptimg2 = sptimg2 / sum(sptimg2(:)) * Isig;
        bgimg   = ones(n1, n3) * Ibg / n1 / n3;

        sptimg4(:,:,n)         = poissrnd(sptimg2) + poissrnd(bgimg) + (100 + 3*randn(n1, n3));
        First_four_rows(:,:,n) = sptimg4(1:4,:,n);
    end
    sptimg4 = uint16(sptimg4) - uint16(repmat(mean(First_four_rows, 3), 4, 1));

    % --- sweep Width x Center ---
    sptn = squeeze(mean(double(sptimg4(7:10,:,:)), 1));   
    for k = 1:num_W
        for j = 1:num_C
            wl_analysis = (Center_Wavelengths(j) - Width_Vals(k)/2) : 1 : (Center_Wavelengths(j) + Width_Vals(k)/2);
            [Std_3D(k,i,j), Accuracy_3D(k,i,j)] = ...
                window_centroid(wl_fixed, sptn, wl_analysis, GT_Centroid_Fixed);
        end
    end

    run_count = run_count + 1;
    waitbar(run_count/total_runs, f1, sprintf('Psig %d (%d/%d)', ...
        Psig_File_Vals(i), run_count, total_runs));
end
close(f1);
PercentError_3D = (Accuracy_3D / GT_Centroid_Fixed) * 100;
disp('Part 1 complete.');

plot_part1_heatmaps(Std_3D, PercentError_3D, Center_Wavelengths, Width_Vals, ...
    Psig_File_Vals, GT_Centroid_Fixed, num_W, font_title, font_label, font_tick, font_cb);

%==========================================================================
% PART 2 : SPECTRAL HETEROGENEITY 
%==========================================================================
if ENABLE_HETERO
spt_raw = Alexa_Fluor_647(:,3); spt_raw = spt_raw(:);
wl      = Alexa_Fluor_647(:,1); wl      = wl(:);

Std_4D      = zeros(num_S, num_W, num_P, num_C);
Accuracy_4D = zeros(num_S, num_W, num_P, num_C);

disp('Part 2 : one realization per (Sigma, Psig), sweep all windows...');
total_runs = num_S * num_P;  run_count = 0;
f2 = waitbar(0, 'Part 2 ...');

for si = 1:num_S
    current_sigma = Sigma_Shift_Vals(si);
    for i = 1:num_P
        Isig = Psig_File_Vals(i) * eta;

        sptimg4         = zeros(16, n3, numSpectra);
        First_four_rows = zeros(4,  n3, numSpectra);
        GT_per_molecule = zeros(1, numSpectra);
        for n = 1:numSpectra
            if current_sigma > 0
                spt_hetero = interp1(wl + randn*current_sigma, spt_raw, wl, 'linear', 0);
                spt_hetero = max(spt_hetero, 0);
            else
                spt_hetero = spt_raw;
            end
            GT_per_molecule(n) = sum(wl .* spt_hetero) / sum(spt_hetero);

            sptimg = conv2(psf, spt_hetero');
            sptimg(:, end-(c1-2):end) = []; sptimg(:, 1:c1) = [];
            sptimg2 = interp1(wl, sptimg', wl_fixed)';
            sptimg2 = sptimg2 / sum(sptimg2(:)) * Isig;
            bgimg   = ones(n1, n3) * Ibg / n1 / n3;

            sptimg4(:,:,n)         = poissrnd(sptimg2) + poissrnd(bgimg) + (100 + 3*randn(n1, n3));
            First_four_rows(:,:,n) = sptimg4(1:4,:,n);
        end
        sptimg4 = uint16(sptimg4) - uint16(repmat(mean(First_four_rows, 3), 4, 1));

        sptn = squeeze(mean(double(sptimg4(7:10,:,:)), 1));
        for k = 1:num_W
            for j = 1:num_C
                wl_analysis = (Center_Wavelengths(j) - Width_Vals(k)/2) : 1 : (Center_Wavelengths(j) + Width_Vals(k)/2);
                [Std_4D(si,k,i,j), Accuracy_4D(si,k,i,j)] = ...
                    window_centroid(wl_fixed, sptn, wl_analysis, GT_per_molecule);
            end
        end

        run_count = run_count + 1;
        waitbar(run_count/total_runs, f2, sprintf('sigma %.1f, Psig %d', ...
            current_sigma, Psig_File_Vals(i)));
    end
end
close(f2);
PercentError_4D = (Accuracy_4D / GT_Centroid_Fixed) * 100;
disp('Part 2 complete.');

plot_part2_figures(Std_4D, PercentError_4D, Sigma_Shift_Vals, Center_Wavelengths, ...
    Width_Vals, Psig_File_Vals, GT_Centroid_Fixed, num_S, num_W, ...
    font_title, font_label, font_tick);
else
    disp('ENABLE_HETERO = false: Part 2 (heterogeneity) skipped.');
end

%==========================================================================
% Functions
%==========================================================================
function [stdC, accC] = window_centroid(wl_fixed, sptn, wl_analysis, GT)
% Vectorized per-molecule centroid over one analysis window.
    if length(wl_fixed) ~= size(sptn,1)
        wl_fixed = wl_fixed(1:size(sptn,1));
    end
    vq = interp1(wl_fixed, sptn, wl_analysis);   
    vq(isnan(vq)) = 0;
    colsum   = sum(vq, 1);
    Centroid = (wl_analysis(:).' * vq) ./ colsum;
    Centroid(colsum <= 0) = NaN;
    stdC = std(Centroid, 'omitnan');
    accC = mean(Centroid - GT, 'omitnan');       
end

function plot_part1_heatmaps(Std_3D, PercentError_3D, Center_Wavelengths, Width_Vals, ...
        Psig_File_Vals, GT_Centroid_Fixed, num_W, font_title, font_label, font_tick, font_cb)
    xtick_wl     = 575:50:775;
    xtick_idx    = arrayfun(@(x) find(Center_Wavelengths==x,1), xtick_wl, 'UniformOutput', false);
    vm           = ~cellfun(@isempty, xtick_idx);
    xtick_idx    = cell2mat(xtick_idx(vm));
    xtick_labels = cellstr(num2str(xtick_wl(vm)'));
    ytick_labels = cellstr(num2str(Width_Vals'));

    Psig_select     = [500, 1000, 2100, 3000, 4600];
    Psig_select_idx = arrayfun(@(p) find(abs(Psig_File_Vals-p)==min(abs(Psig_File_Vals-p)),1), Psig_select);
    Psig_titles     = arrayfun(@(p) sprintf('P_{sig} = %d', p), Psig_select, 'UniformOutput', false);
    ns              = numel(Psig_select);

    [~, c_idx_GT]  = min(abs(Center_Wavelengths - GT_Centroid_Fixed));
    [~, w_idx_75]  = min(abs(Width_Vals - 150));
    [~, w_idx_100] = min(abs(Width_Vals - 200));

    figure('Name','Heatmaps: Std & |% Error|','Position',[50 50 ns*320+100 900]);
    for s = 1:ns
        subplot(2, ns, s);
        imagesc(squeeze(Std_3D(:, Psig_select_idx(s), :)));
        colormap(jet); h = colorbar; ylabel(h, '\sigma_C (nm)', 'FontSize', font_cb);
        hold on;
        plot(c_idx_GT, w_idx_75,'k*','MarkerSize',10,'LineWidth',2);
        plot(c_idx_GT, w_idx_100,'ko','MarkerSize',10,'LineWidth',2); hold off;
        title({'Std Heatmap (Precision)'; Psig_titles{s}}, 'FontSize', font_title);
        xlabel('Window Center (nm)', 'FontSize', font_label);
        if s==1, ylabel('Window Width (nm)', 'FontSize', font_label); end
        set(gca,'XTick',xtick_idx,'XTickLabel',xtick_labels,'FontSize',font_tick);
        set(gca,'YTick',1:num_W,'YTickLabel',ytick_labels,'FontSize',font_tick,'YDir','normal');
        if s==1, legend({'\pm75nm (W=150)','\pm100nm (W=200)'},'Location','northwest','FontSize',font_tick); end

        subplot(2, ns, ns + s);
        imagesc(abs(squeeze(PercentError_3D(:, Psig_select_idx(s), :))));
        colormap(jet); h = colorbar; ylabel(h, '|% Error| (%)', 'FontSize', font_cb);
        hold on;
        plot(c_idx_GT, w_idx_75,'k*','MarkerSize',10,'LineWidth',2);
        plot(c_idx_GT, w_idx_100,'ko','MarkerSize',10,'LineWidth',2); hold off;
        title({'|Percent Error| Heatmap'; Psig_titles{s}}, 'FontSize', font_title);
        xlabel('Window Center (nm)', 'FontSize', font_label);
        if s==1, ylabel('Window Width (nm)', 'FontSize', font_label); end
        set(gca,'XTick',xtick_idx,'XTickLabel',xtick_labels,'FontSize',font_tick);
        set(gca,'YTick',1:num_W,'YTickLabel',ytick_labels,'FontSize',font_tick,'YDir','normal');
        if s==1, legend({'\pm75nm (W=150)','\pm100nm (W=200)'},'Location','northwest','FontSize',font_tick); end
    end
    exportgraphics(gcf, 'Heatmaps_5Psig.pdf', 'ContentType', 'vector');

    fprintf('\n===== Manuscript Range Verification (AF647, all windows) =====\n');
    for s = 1:ns
        sm = squeeze(Std_3D(:, Psig_select_idx(s), :));
        pm = abs(squeeze(PercentError_3D(:, Psig_select_idx(s), :)));
        fprintf('Psig = %-5d | Std: %.2f - %.2f nm | |%%Err|: %.1f - %.1f%%\n', ...
            Psig_select(s), min(sm(:)), max(sm(:)), min(pm(:)), max(pm(:)));
    end
end

function plot_part2_figures(Std_4D, PercentError_4D, Sigma_Shift_Vals, Center_Wavelengths, ...
        Width_Vals, Psig_File_Vals, GT_Centroid_Fixed, num_S, num_W, font_title, font_label, font_tick)
    colors = lines(num_S);
    styles = {'-o','-s','-^','-d'};
    slbl   = arrayfun(@(s) sprintf('\\sigma = %.1f nm', s), Sigma_Shift_Vals, 'UniformOutput', false);

    xtick_wl     = 575:50:775;
    xtick_idx    = arrayfun(@(x) find(Center_Wavelengths==x,1), xtick_wl, 'UniformOutput', false);
    vm           = ~cellfun(@isempty, xtick_idx);
    xtick_idx    = cell2mat(xtick_idx(vm));
    xtick_labels = cellstr(num2str(xtick_wl(vm)'));
    ytick_labels = cellstr(num2str(Width_Vals'));

    [~, c_idx_GT]  = min(abs(Center_Wavelengths - GT_Centroid_Fixed));
    [~, w_idx_150] = min(abs(Width_Vals - 150));
    [~, w_idx_200] = min(abs(Width_Vals - 200));

    fig4 = figure('Name','Heterogeneity Effect','Position',[100 100 1400 1000]);
    wc = [w_idx_150, w_idx_200];
    wl_lbl = {'Width = 150 nm (\pm75nm)','Width = 200 nm (\pm100nm)'};
    for ww = 1:2
        subplot(2,2,(ww-1)*2+1); hold on;
        for si = 1:num_S
            plot(Psig_File_Vals, squeeze(Std_4D(si,wc(ww),:,c_idx_GT)), styles{si}, ...
                'Color', colors(si,:), 'LineWidth',2,'MarkerSize',6,'DisplayName',slbl{si});
        end
        hold off; xlabel('P_{sig} (photons)','FontSize',font_label);
        ylabel('Std \sigma_C (nm)','FontSize',font_label);
        title({'Precision vs P_{sig}'; wl_lbl{ww}},'FontSize',font_title);
        legend('Location','best','FontSize',font_tick);
        set(gca,'XTick',Psig_File_Vals,'FontSize',font_tick); xtickangle(45); grid on;

        subplot(2,2,(ww-1)*2+2); hold on;
        for si = 1:num_S
            plot(Psig_File_Vals, abs(squeeze(PercentError_4D(si,wc(ww),:,c_idx_GT))), styles{si}, ...
                'Color', colors(si,:), 'LineWidth',2,'MarkerSize',6,'DisplayName',slbl{si});
        end
        hold off; xlabel('P_{sig} (photons)','FontSize',font_label);
        ylabel('|Percent Error| (%)','FontSize',font_label);
        title({'|% Error| vs P_{sig}'; wl_lbl{ww}},'FontSize',font_title);
        legend('Location','best','FontSize',font_tick);
        set(gca,'XTick',Psig_File_Vals,'FontSize',font_tick); xtickangle(45); grid on;
    end
    exportgraphics(fig4, 'Figure4_Contrast_Plot.pdf', 'ContentType', 'vector');

    [~, ip] = min(abs(Psig_File_Vals - 2100));
    fig5s = figure('Name','Heterogeneity Heatmaps: Std','Position',[100 100 1600 1200]);
    for si = 1:num_S
        subplot(2,2,si); imagesc(squeeze(Std_4D(si,:,ip,:))); colormap(jet); colorbar;
        title({sprintf('Std Heatmap, \\sigma = %.1f nm', Sigma_Shift_Vals(si)); ...
               sprintf('P_{sig} = %d', Psig_File_Vals(ip))},'FontSize',font_title);
        xlabel('Window Center (nm)','FontSize',font_label); ylabel('Window Width (nm)','FontSize',font_label);
        set(gca,'XTick',xtick_idx,'XTickLabel',xtick_labels,'FontSize',font_tick);
        set(gca,'YTick',1:num_W,'YTickLabel',ytick_labels,'FontSize',font_tick,'YDir','normal');
    end
    exportgraphics(fig5s, 'Figure5_Heatmap_Std.pdf', 'ContentType', 'vector');

    fig5e = figure('Name','Heterogeneity Heatmaps: |% Error|','Position',[100 100 1600 1200]);
    for si = 1:num_S
        subplot(2,2,si); imagesc(abs(squeeze(PercentError_4D(si,:,ip,:)))); colormap(jet); colorbar;
        title({sprintf('|%%Error| Heatmap, \\sigma = %.1f nm', Sigma_Shift_Vals(si)); ...
               sprintf('P_{sig} = %d', Psig_File_Vals(ip))},'FontSize',font_title);
        xlabel('Window Center (nm)','FontSize',font_label); ylabel('Window Width (nm)','FontSize',font_label);
        set(gca,'XTick',xtick_idx,'XTickLabel',xtick_labels,'FontSize',font_tick);
        set(gca,'YTick',1:num_W,'YTickLabel',ytick_labels,'FontSize',font_tick,'YDir','normal');
    end
    exportgraphics(fig5e, 'Figure5_Heatmap_Error.pdf', 'ContentType', 'vector');

    [~, c700] = min(abs(Center_Wavelengths - 700));
    [~, wopt] = min(abs(Width_Vals - 150));
    figure('Name','Sweet Spot vs Heterogeneity','Position',[100 100 1200 500]);
    subplot(1,2,1); hold on;
    for si = 1:num_S
        plot(Psig_File_Vals, squeeze(Std_4D(si,wopt,:,c700)), styles{si}, ...
            'Color', colors(si,:),'LineWidth',2,'MarkerSize',6,'DisplayName',slbl{si});
    end
    yline(1.5,'k--','Target 1.5 nm','LineWidth',1.5,'FontSize',font_tick); hold off;
    xlabel('P_{sig} (photons)','FontSize',font_label); ylabel('Std \sigma_C (nm)','FontSize',font_label);
    title({'Precision vs P_{sig}';'Center = 700 nm, Width = 150 nm'},'FontSize',font_title);
    legend('Location','best','FontSize',font_tick);
    set(gca,'XTick',Psig_File_Vals,'FontSize',font_tick); xtickangle(45); grid on;

    subplot(1,2,2); hold on;
    for si = 1:num_S
        plot(Psig_File_Vals, abs(squeeze(PercentError_4D(si,wopt,:,c700))), styles{si}, ...
            'Color', colors(si,:),'LineWidth',2,'MarkerSize',6,'DisplayName',slbl{si});
    end
    hold off; xlabel('P_{sig} (photons)','FontSize',font_label);
    ylabel('|Percent Error| (%)','FontSize',font_label);
    title({'|% Error| vs P_{sig}';'Center = 700 nm, Width = 150 nm'},'FontSize',font_title);
    legend('Location','best','FontSize',font_tick);
    set(gca,'XTick',Psig_File_Vals,'FontSize',font_tick); xtickangle(45); grid on;
end
