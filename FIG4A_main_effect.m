clear;clc;

%% 1. load data
load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat')
Demo_HCs = readtable('/host/verges/tank/data/mo/BIDS_CERMEP-MRXFDG/BIDS_CERMEP-MRXFDG_demo.csv');
Demo_PTs = readtable('/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/dis_BTTH1_BFTH.csv');

Age_HCs = Demo_HCs.Age; Age_PTs = Demo_PTs.AgeSurg; Age_all = [Age_HCs; Age_PTs]; Prog_PTs = Demo_PTs.x1YProgI;
Sex_HCs = Demo_HCs.Sex; Sex_PTs = Demo_PTs.Sex;     Sex_all = [Sex_HCs; Sex_PTs]; Sex_all_num = double(strcmp(Sex_all,'M'));
Dur_PTs = Demo_PTs.Duration;
SZFreq_PTs = Demo_PTs.SZFreq;   SZFreq_PTs_num = zeros(size(SZFreq_PTs));
SZFreq_PTs_num(strcmp(SZFreq_PTs, 'Daily')) = 1;   SZFreq_PTs_num(strcmp(SZFreq_PTs, 'Weekly')) = 2;
SZFreq_PTs_num(strcmp(SZFreq_PTs, 'Monthly')) = 3; SZFreq_PTs_num(strcmp(SZFreq_PTs, 'Yearly')) = 4;

% % I already run and save the data of Mahalanubis Distance
% %% 2. Compute Mahalanubis Distance
% HippoThick_HCs  = Voxel_HCs_subfeats_hmz_flip_asym_wPDB.HippoThick;
% HippoCurv_HCs   = Voxel_HCs_subfeats_hmz_flip_asym_wPDB.HippoCurv;
% HippoFLAIR_HCs  = Voxel_HCs_subfeats_hmz_flip_asym_wPDB.HippoFLAIR;
% HippoMyelin_HCs = Voxel_HCs_subfeats_hmz_flip_asym_wPDB.HippoMyelin;
% HippoPET_HCs    = Voxel_HCs_subfeats_hmz_flip_asym_wPDB.HippoPET;
% 
% HippoThick_PTs  = Voxel_PTs_subfeats_hmz_flip_asym_wPDB.HippoThick;
% HippoCurv_PTs   = Voxel_PTs_subfeats_hmz_flip_asym_wPDB.HippoCurv;
% HippoFLAIR_PTs  = Voxel_PTs_subfeats_hmz_flip_asym_wPDB.HippoFLAIR;
% HippoMyelin_PTs = Voxel_PTs_subfeats_hmz_flip_asym_wPDB.HippoMyelin;
% HippoPET_PTs    = Voxel_PTs_subfeats_hmz_flip_asym_wPDB.HippoPET;
% 
% MidPET_HCs      = Voxel_HCs_neofeats_hmz_flip_asym_wPDB.MidPET;
% MidPET_PTs      = Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET;
% 
% % compute Mahalanobis Distance
% wPDB_COV_HCs = [Age_HCs(:), Sex_all_num(1:length(Demo_HCs.ID))];
% wPDB_COV_PTs = [Age_PTs(:), Sex_all_num(length(Demo_HCs.ID)+1:end)];
% 
% % Combine all features
% X_HCs = [HippoThick_HCs, HippoCurv_HCs, HippoFLAIR_HCs, HippoMyelin_HCs, HippoPET_HCs];
% X_PTs = [HippoThick_PTs, HippoCurv_PTs, HippoFLAIR_PTs, HippoMyelin_PTs, HippoPET_PTs];
% 
% % Get dimensions
% n_regions = size(HippoThick_HCs, 2);
% n_features = 5; % Number of feature types
% n_HCs = size(X_HCs, 1);
% n_PTs = size(X_PTs, 1);
% 
% % Preallocate
% Mahal_HCs = zeros(n_HCs, n_regions);
% Mahal_PTs = zeros(n_PTs, n_regions);
% 
% % parpool
% 
% parfor region = 1:n_regions
%     % Extract features for this region across all feature types
%     region_indices = region:n_regions:(n_features*n_regions);
%     X_HCs_region = X_HCs(:, region_indices);
%     X_PTs_region = X_PTs(:, region_indices);
% 
%     % Control mean and covariance
%     mu = mean(X_HCs_region, 1);
%     Sigma = cov(X_HCs_region);
%     invSigma = pinv(Sigma);  % safer than inv()
% 
%     % Mahalanobis distance for Controls
%     diff_HCs = X_HCs_region - mu;
%     Mahal_HCs_temp = sqrt(sum((diff_HCs * invSigma) .* diff_HCs, 2));
% 
%     % Mahalanobis distance for Patients
%     diff_PTs = X_PTs_region - mu;
%     Mahal_PTs_temp = sqrt(sum((diff_PTs * invSigma) .* diff_PTs, 2));
% 
%     % Store results (this indexing is allowed in parfor)
%     Mahal_HCs(:, region) = Mahal_HCs_temp;
%     Mahal_PTs(:, region) = Mahal_PTs_temp;
% end
% 
% % W-score calculation
% Mahal_wPDB_HCs = zeros(n_HCs, n_regions);
% Mahal_wPDB_PTs = zeros(n_PTs, n_regions);
% 
% for region = 1:n_regions
%     Vals_HCs = Mahal_HCs(:, region);
%     Vals_PTs = Mahal_PTs(:, region);
% 
%     % W-score using the same covariates as before
%     [Mahal_wPDB_HCs(:, region), Mahal_wPDB_PTs(:, region)] = ...
%         wscore_all(Vals_HCs, wPDB_COV_HCs, Vals_PTs, wPDB_COV_PTs);
% end
% 
% save('/host/verges/tank/data/mo/demo/31Voxel_Mahal_wPDB', ...
%     'Mahal_wPDB_HCs', 'Mahal_wPDB_PTs')

%% 3. load_data
load_color
load_surface
load('/host/verges/tank/data/mo/demo/31Voxel_Mahal_wPDB.mat')
load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat')

%% 4. plot the folded and unfolded surface of all features
hippoSurf_lhFile       = gifti('/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/derivatives/hippunfold_v1.4.1/hippunfold/sub-MNI152NLin2009cAsym/surf/sub-MNI152NLin2009cAsym_hemi-L_space-T1w_den-2mm_label-hipp_midthickness.surf.gii');
hippoUnfoldSurf_lhFile = gifti('/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/derivatives/hippunfold_v1.4.1/hippunfold/sub-MNI152NLin2009cAsym/surf/sub-MNI152NLin2009cAsym_hemi-L_space-unfold_den-2mm_label-hipp_midthickness.surf.gii');

hippoSurf_rhFile       = gifti('/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/derivatives/hippunfold_v1.4.1/hippunfold/sub-MNI152NLin2009cAsym/surf/sub-MNI152NLin2009cAsym_hemi-R_space-T1w_den-2mm_label-hipp_midthickness.surf.gii');
hippoUnfoldSurf_rhFile = gifti('/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/derivatives/hippunfold_v1.4.1/hippunfold/sub-MNI152NLin2009cAsym/surf/sub-MNI152NLin2009cAsym_hemi-R_space-unfold_den-2mm_label-hipp_midthickness.surf.gii');

hippoSurf       = struct();
hippoSurf.coord = [hippoSurf_lhFile.vertices', hippoSurf_rhFile.vertices'];
hippoSurf.tri   = [hippoSurf_lhFile.faces;     hippoSurf_rhFile.faces + 419];

hippoUnfoldSurf = struct();
hippoUnfoldSurf.coord = [hippoUnfoldSurf_lhFile.vertices',  hippoUnfoldSurf_rhFile.vertices'];
hippoUnfoldSurf.tri   = [hippoUnfoldSurf_lhFile.faces;      hippoUnfoldSurf_rhFile.faces + 419];

%  4.1 bkg data
figure
bkgData = ones(1,838);
SurfStatViewDataHippo_Unfold(bkgData, hippoSurf, hippoUnfoldSurf, 'bkg', 'white')
colormap([0 0 0])

%  4.2 Mahalanobis distance
figure
MahalData = double(mean(Mahal_wPDB_PTs,1)');
outlierIdx = isoutlier(MahalData);
MahalData_clean = MahalData;
inIdx  = ~outlierIdx & ~isnan(MahalData);
fillIdx = outlierIdx | isnan(MahalData);
V = double(hippoSurf.coord');
F = scatteredInterpolant(V(inIdx,1), V(inIdx,2), V(inIdx,3), ...
                          double(MahalData(inIdx)), ...
                          'natural','nearest');
MahalData_clean(fillIdx) = F(V(fillIdx,1), ...
                             V(fillIdx,2), ...
                             V(fillIdx,3) );
SurfStatViewDataHippo_Unfold(MahalData_clean, hippoSurf, hippoUnfoldSurf, 'Mahalanobis', 'white')
colormap(red)

%  4.3 Each feature
feature_names = {'HippoThick', 'HippoCurv', 'HippoFLAIR', 'HippoMyelin', 'HippoPET'};
cmaps_color = {cmap_purple, cmap_blue, cmap_green, cmap_orange, cmap_yellow};
for i = 1:numel(feature_names)
    figure
    featData = double(mean(Voxel_PTs_subfeats_hmz_flip_asym_wPDB.(feature_names{i}),1)');  % N×1
    out = isoutlier(featData);
    featData_clean = featData;
    inIdx   = ~out;          % or: ~out & ~isnan(x)
    fillIdx = out;           % or: out | isnan(x)
    F = scatteredInterpolant(V(inIdx,1), V(inIdx,2), V(inIdx,3), ...
                             featData(inIdx), 'natural','nearest');
    featData_clean(fillIdx) = F(V(fillIdx,1), V(fillIdx,2), V(fillIdx,3));
    SurfStatViewDataHippo_Unfold(featData_clean, hippoSurf, hippoUnfoldSurf, feature_names{i}, 'white');
    colormap(cmaps_color{i});
    clim([0 100])
end

%% 5. plot the hippoFeat on neocortex
Age_Cov = term(Age_PTs); Sex_Cov    = term(Sex_PTs);
Dur_Cov = term(Dur_PTs); SZFreq_Cov = term(SZFreq_PTs_num);

%  5.1 Mahal
medialIdx = find(economo_parc(1:end/2,1) ~= 0);
PET_lh = Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET(:, 1:end/2);
PET_lh(:, setdiff(1:size(PET_lh,2), medialIdx)) = 0;

Madal_wscore_PTs_mean = mean(Mahal_wPDB_PTs(:,1:end/2),2);   Madal_Var = term(Madal_wscore_PTs_mean);

Model = 1 + Madal_Var + Age_Cov + Sex_Cov + Dur_Cov + SZFreq_Cov;
slm_Madal = SurfStatLinMod(PET_lh, Model);
slm_Madal_Pos = SurfStatT(slm_Madal, Madal_wscore_PTs_mean);
slm_Madal_Neg = SurfStatT(slm_Madal, -Madal_wscore_PTs_mean);

p_Pos = 1-tcdf(slm_Madal_Pos.t,slm_Madal_Pos.df);
p_Neg = 1-tcdf(slm_Madal_Neg.t,slm_Madal_Neg.df);
[h_Pos, crit_p_Pos, adj_p_Pos] = fdr_bh(p_Pos,0.05);
[h_Neg, crit_p_Neg, adj_p_Neg] = fdr_bh(p_Neg,0.05);

% raw cohen's D    % sig_mask = (p_Pos < 0.05) | (p_Neg < 0.05);
D_Madal = slm_Madal_Neg.ef;

sigMask_Mahal_uncorr = (p_Neg < 0.05);
D_Madal_uncorr = D_Madal .* sigMask_Mahal_uncorr;
plot_hemispheres(D_Madal_uncorr', {surf_lh}); colormap(cmap_red); 

% after FDR
sigMask_Mahal_FDR = (h_Neg == 1);
D_Madal_FDR = D_Madal .* sigMask_Mahal_FDR;
plot_hemispheres(D_Madal_FDR', {surf_lh}); colormap(flipud(cmap_red)); 

%  5.2 plot individual neo features for the negative results
for i = 1:length(feature_names)
    current_feature = Voxel_PTs_subfeats_hmz_flip_asym_wPDB.(feature_names{i});
    Hippo_feats_temp = mean(current_feature,2);
    Hippo_feats_temp_Var = term(Hippo_feats_temp);
    
    Model = 1 + Hippo_feats_temp_Var + Age_Cov + Sex_Cov + Dur_Cov + SZFreq_Cov;
    slm_temp = SurfStatLinMod(PET_lh, Model);
    slm_temp_Pos = SurfStatT(slm_temp, Hippo_feats_temp);
    slm_temp_Neg = SurfStatT(slm_temp, -Hippo_feats_temp);
    
    p_temp_Pos = 1 - tcdf(slm_temp_Pos.t, slm_temp_Pos.df);
    p_temp_Neg = 1 - tcdf(slm_temp_Neg.t, slm_temp_Neg.df);
      
    [h_Pos, crit_p_Pos, adj_p_Pos] = fdr_bh(p_temp_Pos, 0.05);
    [h_Neg, crit_p_Neg, adj_p_Neg] = fdr_bh(p_temp_Neg, 0.05);
    
    D_temp = slm_temp_Neg.ef;

    sigMask_uncorr = (p_temp_Neg < 0.05);
    D_temp_uncorr(i,:) = D_temp .* sigMask_uncorr;
    plot_hemispheres(D_temp_uncorr(i,:)', {surf_lh}, 'labeltext', feature_names{1,i}); colormap(cmap_blue);
    max(D_temp_uncorr(i,:))

    h_Neg_temp = double(h_Neg);
    h_Neg_temp(h_Neg_temp == 0) = 1;
    plot_hemispheres(h_Neg_temp', {surf_lh}, 'labeltext', feature_names{1,i}); colormap(gray);
    set(findall(gcf, 'type', 'axes'), 'CLim', [0 1]);

    D_hippoFeats(i,:) = D_temp;
end

%  5.2 plot individual neo features for the positive results
for i = 1:length(feature_names)
    current_feature = Voxel_PTs_subfeats_hmz_flip_asym_wPDB.(feature_names{i});
    Hippo_feats_temp = mean(current_feature,2);
    Hippo_feats_temp_Var = term(Hippo_feats_temp);
    
    Model = 1 + Hippo_feats_temp_Var + Age_Cov + Sex_Cov + Dur_Cov + SZFreq_Cov;
    slm_temp = SurfStatLinMod(PET_lh, Model);
    slm_temp_Pos = SurfStatT(slm_temp, Hippo_feats_temp);
    slm_temp_Neg = SurfStatT(slm_temp, -Hippo_feats_temp);
    
    p_temp_Pos = 1 - tcdf(slm_temp_Pos.t, slm_temp_Pos.df);
    p_temp_Neg = 1 - tcdf(slm_temp_Neg.t, slm_temp_Neg.df);
      
    [h_Pos, crit_p_Pos, adj_p_Pos] = fdr_bh(p_temp_Pos, 0.001);
    [h_Neg, crit_p_Neg, adj_p_Neg] = fdr_bh(p_temp_Neg, 0.001);
    
    D_temp = slm_temp_Pos.ef;

    sigMask_uncorr = (p_temp_Pos < 0.05);
    D_temp_uncorr(i,:) = D_temp .* sigMask_uncorr;
    plot_hemispheres(D_temp_uncorr(i,:)', {surf_lh}, 'labeltext', feature_names{1,i}); colormap(cmap_blue);
    
    h_Pos_temp = double(h_Pos);
    h_Pos_temp(h_Pos_temp == 0) = -1;
    plot_hemispheres(h_Pos_temp', {surf_lh}, 'labeltext', feature_names{1,i}); colormap(gray);
    set(findall(gcf, 'type', 'axes'), 'CLim', [-1 0]);

    D_hippoFeats(i,:) = D_temp;
end

% save('/host/verges/tank/data/mo/demo/42Voxel_Dvalue', ...
%     'D_hippoFeats', 'D_Madal')

%% plot the distribution of HippoFeats (Boxplot) for PTs
rng(42);

var_names = {'HippoThick','HippoCurv','HippoFLAIR','HippoMyelin','HippoPET','Mahalanobis'};
cmaps     = {cmap_purple, cmap_blue, cmap_green, cmap_orange, cmap_yellow, cmap_red};
n_scatter = 500;   % subsampled scatter points
bw        = 0.12;  % narrower box
x_pos     = 1;

load_paper

for vi = 1:6
    ax = subplot(10, 7, 7*vi);
    hold(ax, 'on');

    % ── Load data ──────────────────────────────────────────────────────
    if vi < 6
        featData_temp = double(mean(Voxel_PTs_subfeats_hmz_flip_asym_wPDB.(var_names{vi}), 1)');
        featData = featData_temp(1:end/2,1);
    else
        featData = double(mean(Mahal_wPDB_PTs(:, 1:end/2), 2));
    end

    % ── 1) Remove outliers (GESD method) ──────────────────────────────
    keep     = ~isoutlier(featData, 'gesd');
    featData = featData(keep);

    cmap_v  = cmaps{vi};
    color_g = cmap_v(180, :);
    solid_c = [0 0 0];   % black line color

    % ── Box stats ──────────────────────────────────────────────────────
    q1    = prctile(featData, 25);
    q3    = prctile(featData, 75);
    mn    = mean(featData);       % mean line
    iqr_v = q3 - q1;
    w_lo  = max(featData(featData >= q1 - 1.5*iqr_v));
    w_hi  = min(featData(featData <= q3 + 1.5*iqr_v));

    % ── Scatter BEHIND box (jittered, subsampled) ────────────────────
    n_plot = min(n_scatter, numel(featData));
    s_idx  = randperm(numel(featData), n_plot);
    jitter = (rand(n_plot, 1) - 0.5) * bw * 1.8;
    scatter(ax, x_pos + jitter, featData(s_idx), 4, ...
            color_g, 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');

    % ── Box fill ────────────────────────────────────────────────────────
    fill(ax, [x_pos-bw, x_pos+bw, x_pos+bw, x_pos-bw], ...
             [q1, q1, q3, q3], color_g, ...
             'EdgeColor', solid_c, 'LineWidth', 0.5, 'FaceAlpha', 0.75);

    % ── Mean line ───────────────────────────────────────────────────────
    plot(ax, [x_pos-bw, x_pos+bw], [mn, mn], '-', ...
         'Color', solid_c, 'LineWidth', 2);

    % ── Whiskers ────────────────────────────────────────────────────────
    plot(ax, [x_pos, x_pos], [w_lo, q1], '-', 'Color', solid_c, 'LineWidth', 0.5);
    plot(ax, [x_pos, x_pos], [q3, w_hi], '-', 'Color', solid_c, 'LineWidth', 0.5);

    % ── Whisker caps ────────────────────────────────────────────────────
    plot(ax, [x_pos-bw*0.5, x_pos+bw*0.5], [w_lo, w_lo], '-', 'Color', solid_c, 'LineWidth', 0.5);
    plot(ax, [x_pos-bw*0.5, x_pos+bw*0.5], [w_hi, w_hi], '-', 'Color', solid_c, 'LineWidth', 0.5);

    % ── Auto ylim with padding top and bottom ───────────────────────────
    pad      = (max(featData) - min(featData)) * 0.15;
    pad      = max(pad, 0.2);
    y_lo_lim = min(featData) - pad;
    y_hi_lim = max(featData) + pad;

    % ── Formatting ──────────────────────────────────────────────────────
    box(ax, 'off');
    xlim(ax, [0.6, 1.4]);
    ylim(ax, [y_lo_lim, y_hi_lim]);
    xticks(ax, []);
    ylabel(ax, var_names{vi}, 'Interpreter', 'none', 'FontSize', 8);

    % Auto y-ticks: low, mid, high rounded to 1 decimal
    y_mid = round((y_lo_lim + y_hi_lim) / 2, 1);
    set(ax, 'YTick', [round(y_lo_lim,1), y_mid, round(y_hi_lim,1)], 'FontSize', 7);

    % ── Slim down subplot width via Position ────────────────────────────
    pos = get(ax, 'Position');
    set(ax, 'Position', [pos(1)+pos(3)*0.15, pos(2), pos(3)*0.6, pos(4)]);
end