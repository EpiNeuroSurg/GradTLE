clear;clc;

%% 1. Load data
load_demo
load_surface
load_color
load_brh

%% 2. Surgical metrics surface representation
load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat')
PATH_PTs = '/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/derivatives/';

for i = 153
    %  3.1 surg 
    File_PTs_LH_surgMask = strcat(PATH_PTs,'surgMask/',Demo_disPTs.ID{i},'/ses-post/',Demo_disPTs.ID{i},'_ses-post_hemi-L_space-nativepro_surf-fsLR-32k_label-cavity.func.gii');
    File_PTs_RH_surgMask = strcat(PATH_PTs,'surgMask/',Demo_disPTs.ID{i},'/ses-post/',Demo_disPTs.ID{i},'_ses-post_hemi-R_space-nativepro_surf-fsLR-32k_label-cavity.func.gii');
    assume_surg_raw = [gifti(File_PTs_LH_surgMask).cdata; gifti(File_PTs_RH_surgMask).cdata];
    
    if strcmp(Demo_disPTs.Side{i}, 'R')
    assume_surg_temp = [assume_surg_raw(end/2+1:end, 1); assume_surg_raw(1:end/2, 1)];
    else
    assume_surg_temp = assume_surg_raw;
    end
    
    in_mask = ismember(economo_parc, allTempLbls);  
    assume_surg_temp(~in_mask) = 0;  
    assume_surg = assume_surg_temp(1:end/2, 1); 
    % plot_hemispheres(assume_surg, surf_lh, 'views', 'lm', 'labeltext', {Demo_disPTs.ID(i)});     colormap(white2black);

    % 3.2 indi PET
    indi_MidPET_hmz_flip_asym_wPDB = Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET(i,1:32492)';
    indi_MidPET_hmz_flip_asym_wPDB(economo_parc(1:end/2,1) == 0) = max(indi_MidPET_hmz_flip_asym_wPDB);
    indi_MidPET_hmz_flip_asym_wPDB(~in_mask) = 0;
    % plot_hemispheres(indi_MidPET_hmz_flip_asym_wPDB(1:end/2, 1), surf_lh, 'views', 'lm', 'labeltext', {Demo_disPTs.ID(i)}); colormap(flipud(blue));

    assume_surg2 = assume_surg;
    assume_surg2(assume_surg2 == 0) = NaN;
    assume_surg2(1,1) = 0;
    % plot_hemispheres(assume_surg2, surf_lh, 'views', 'lm', 'labeltext', {Demo_disPTs.ID(i)}); colormap(white2black);
        
    % 3.3 w-score PET
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal = indi_MidPET_hmz_flip_asym_wPDB;
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal(~(in_mask(1:end/2, 1))) = 0;
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal(indi_MidPET_hmz_flip_asym_wPDB_ThrVal > -1.96) = 0;
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal(indi_MidPET_hmz_flip_asym_wPDB_ThrVal <= -1.96) = 0.8;

    % plot_hemispheres(indi_MidPET_hmz_flip_asym_wPDB_ThrVal(1:end/2, 1), surf_lh, 'views', 'lm', 'labeltext', {Demo_disPTs.ID{i}}); colormap(blue); set(findall(gcf,'type','axes'),'CLim',[0 1]);
end

%% 3. Surgical metrics correlation representation
levels = 0:contour_step_mm:max_dist_mm;
nL     = numel(levels);
half   = numel(dist_roi)/2;

mean_all = nan(nL, 6);
ci_all   = nan(nL, 6);
n_all    = zeros(nL, 6);

vars      = {assume_surg + 3.2, indi_MidPET_hmz_flip_asym_wPDB_ThrVal(1:end/2, 1) + 3.9, indi_MidPET_hmz_flip_asym_wPDB(1:end/2, 1)};
var_names = {'assume_surg', 'indi_ThrVal', 'indi_raw'};
cmaps     = {cmap_mean, cmap_G1, cmap_blue};

% Compute band stats for all 3 vars
for v = 3
    val = vars{v};
    for k = 1:nL
        band = is_TL(1:half) & abs(dist_roi(1:half) - levels(k)) <= contour_tol_mm;
        band = band & ~isnan(val(1:half));
        n_all(k,v) = nnz(band);
        if n_all(k,v) > 1
            s = val(band);
            mean_all(k,v) = mean(s, 'omitnan');
            ci_all(k,v)   = 1.96 * std(s, 'omitnan') / sqrt(n_all(k,v));
        end
    end
end

load_paper
n_colors = 192;
norm_idx = @(n) max(1, min(n_colors, round((n - min(n)) / max(max(n)-min(n), 1) * (n_colors-1) + 1)));

ax1 = subplot(10, 7, [23, 24, 25, 30, 31, 32]);
hold(ax1, 'on');

%% --- v = 1, 2: Binary vars → scatter only where value == 1 ---
for v = 1:2
    cmap_v = cmaps{v}(65:256, :);
    clr    = cmap_v(end, :);
    val    = vars{v};
    d      = dist_roi(1:half);

    % Only keep is_TL points where value == 1 (exclude 0s)
    mask = is_TL(1:half) & val(1:half) > 4 & ~isnan(val(1:half));
    x_sc = d(mask);
    y_sc = val(mask);

    % Small jitter so overlapping dots are visible
    jitter = (rand(size(y_sc)) - 0.5) * 0.04;
    scatter(ax1, x_sc, y_sc + jitter, 8, ...
        'MarkerFaceColor', clr, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.35);

    % Mean proportion line (trend of how many 1s per band)
    idx_v = norm_idx(n_all(:,v));
    valid = ~isnan(mean_all(:,v));
    x_m   = levels(valid)';
    m_m   = mean_all(valid, v);
    idx_s = idx_v(valid);
    for i = 1:numel(x_m)-1
        plot(ax1, x_m(i:i+1), m_m(i:i+1), '-', ...
            'Color', cmap_v(idx_s(i),:), 'LineWidth', 1.5);
    end
end

%% --- v = 3: Continuous var → mean ± 95% CI ribbon ---
v      = 3;
cmap_v = cmaps{v}(65:256, :);
idx_v  = norm_idx(n_all(:,v));
valid  = ~isnan(mean_all(:,v)) & ~isnan(ci_all(:,v));
x      = levels(valid)';
m      = mean_all(valid, v);
ci     = ci_all(valid, v);
idx_s  = idx_v(valid);
clr    = cmap_v(end, :);

fill(ax1, [x; flipud(x)], [m+ci; flipud(m-ci)], clr, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);

for i = 1:numel(x)-1
    plot(ax1, x(i:i+1), m(i:i+1), '-', ...
        'Color', cmap_v(idx_s(i),:), 'LineWidth', 2);
end

box(ax1, 'off');
xlabel(ax1, 'Geodesic distance (mm)');
ylabel(ax1, 'NormPETAsymIdx');
xlim(ax1, [-5 130]);

%% Representation of surface-based data
maskTemp = double(ismember(economo_parc, pureTempLbls));
for i = 65
    % indi PET
    indi_MidPET_hmz_flip_asym_wPDB = Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET(i,1:32492)';
    % plot_hemispheres(indi_MidPET_hmz_flip_asym_wPDB, surf_lh, 'views', 'lm'); colormap(flipud(blue))

    % w-score PETmap
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal = indi_MidPET_hmz_flip_asym_wPDB;
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal(~(maskTemp(1:end/2, 1))) = 0;
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal(indi_MidPET_hmz_flip_asym_wPDB_ThrVal > -1.96) = 0;
    % plot_hemispheres(indi_MidPET_hmz_flip_asym_wPDB_ThrVal, surf_lh, 'views', 'lm'); colormap(flipud(blue)); 

    % surg 
    File_PTs_LH_surgMask = strcat(PATH_PTs,'surgMask/',Demo_disPTs.ID{i},'/ses-post/',Demo_disPTs.ID{i},'_ses-post_hemi-L_space-nativepro_surf-fsLR-32k_label-cavity.func.gii');
    File_PTs_RH_surgMask = strcat(PATH_PTs,'surgMask/',Demo_disPTs.ID{i},'/ses-post/',Demo_disPTs.ID{i},'_ses-post_hemi-R_space-nativepro_surf-fsLR-32k_label-cavity.func.gii');
    assume_surg_raw = [gifti(File_PTs_LH_surgMask).cdata; gifti(File_PTs_RH_surgMask).cdata];

    if strcmp(Demo_disPTs.Side{i}, 'R')
        assume_surg_temp = [assume_surg_raw(end/2+1:end, 1); assume_surg_raw(1:end/2, 1)];
    else
        assume_surg_temp = assume_surg_raw;
    end
    
    in_mask = ismember(economo_parc, allTempLbls);  
    assume_surg_temp(~in_mask) = 0;               
    assume_surg = assume_surg_temp(1:end/2, 1);
    
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal(indi_MidPET_hmz_flip_asym_wPDB_ThrVal ~= 0) = 1;
    % plot_hemispheres([indi_MidPET_hmz_flip_asym_wPDB, indi_MidPET_hmz_flip_asym_wPDB_ThrVal, assume_surg], very_inflated_surf_lh, 'views', 'lm'); colormap(flipud(blue)); 
    plot_hemispheres(indi_MidPET_hmz_flip_asym_wPDB_ThrVal, surf_lh, 'views', 'lm'); colormap([1 1 1; red_mean]);
    plot_hemispheres(assume_surg, surf_lh, 'views', 'lm'); colormap(cmap_purple); 
    plot_hemispheres(assume_surg, surf_lh, 'views', 'lm'); colormap(white2black);
    plot_hemispheres(indi_MidPET_hmz_flip_asym_wPDB, surf_lh, 'views', 'lm'); colormap(flipud(blue)); 

end