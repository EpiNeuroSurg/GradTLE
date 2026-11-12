clear;clc;

%% 1. load data, surface, color
load_demo
load_surface
load_color

%% 2.1 Parcel-wise Cohen's D and ROI plot; Economo
load('/host/verges/tank/data/mo/demo/02Parcel_feats_hmz_flip_asym_wPDB.mat')

allTempLbls = [periTempLbls, mesTempLbls, latTempLbls];
numRegions = length(allTempLbls);
% Two groups only
numGroups = 2;
% Prepare data
% Plot_RegDiff_data = cell(numRegions, numGroups);
for i = 1:numRegions
    Plot_RegDiff_data{i,1} = Parcel_PTs_feats_hmz_flip_asym_wPDB.MidPET(:, allTempLbls(i));
    Plot_RegDiff_data{i,2} = Parcel_HCs_feats_hmz_flip_asym_wPDB.MidPET(:, allTempLbls(i));
end

% figure;
% set(gcf,'Position',[100 100 1600 800]);
% colors = lines(numGroups);
stats = struct;
for f = 1:numRegions
    % subplot(3,5,f); hold on;   
    % Extract data
    group1 = Plot_RegDiff_data{f,1};
    group2 = Plot_RegDiff_data{f,2};
    allData = [group1; group2];
    groupIdx = [ones(size(group1)); 2*ones(size(group2))];
    % Statistics: t-test
    [~,p] = ttest2(group1, group2);
    stats(f).ttest_p = p;
    % Effect size: Cohen's d
    mean1 = mean(group1);
    mean2 = mean(group2);
    sd1 = std(group1);
    sd2 = std(group2);
    pooledSD = sqrt(((length(group1)-1)*sd1^2 + (length(group2)-1)*sd2^2) ...
                    / (length(group1)+length(group2)-2));
    d = (mean1 - mean2) / pooledSD;
    stats(f).cohen_d = d;
    % % Boxplot
    % boxplot(allData, groupIdx, 'Colors','k', 'Symbol','');
    % % Scatter (jitter)
    % x1 = 1 + 0.1*(rand(size(group1))-0.5);
    % x2 = 2 + 0.1*(rand(size(group2))-0.5);
    % scatter(x1, group1, 35, colors(1,:), 'filled', 'MarkerFaceAlpha',0.6);
    % scatter(x2, group2, 35, colors(2,:), 'filled', 'MarkerFaceAlpha',0.6);
    % % Plot formatting
    % xticks([1 2]);
    % xticklabels({'PTs','HCs'});
    % ylabel('Value');
    % grid off;
end
Parcel_D = [stats.cohen_d];
[minParcel_D, idxLowestD] = min(Parcel_D); display([minParcel_D, allTempLbls(1,idxLowestD)]);

%% 2.2 Parcel-wise Cohen's D and ROI plot; Julich
load('/host/verges/tank/data/mo/demo/02Parcel2_feats_hmz_flip_asym_wPDB.mat')
% allTempLbls = [periTempLbls, mesTempLbls, latTempLbls];
numRegions = 148;

% Prepare data
for i = 1:numRegions
    Plot_RegDiff_data{i,1} = Parcel2_PTs_feats_hmz_flip_asym_wPDB.MidPET(:, i);
    Plot_RegDiff_data{i,2} = Parcel2_HCs_feats_hmz_flip_asym_wPDB.MidPET(:, i);
end

% figure;
% set(gcf,'Position',[100 100 1600 800]);
% colors = lines(numGroups);
stats2 = struct;
for f = 1:numRegions
    % subplot(3,5,f); hold on;   
    % Extract data
    group1 = Plot_RegDiff_data{f,1};
    group2 = Plot_RegDiff_data{f,2};
    allData = [group1; group2];
    groupIdx = [ones(size(group1)); 2*ones(size(group2))];
    % Statistics: t-test
    [~,p] = ttest2(group1, group2);
    stats2(f).ttest_p = p;
    % Effect size: Cohen's d
    mean1 = mean(group1);
    mean2 = mean(group2);
    sd1 = std(group1);
    sd2 = std(group2);
    pooledSD = sqrt(((length(group1)-1)*sd1^2 + (length(group2)-1)*sd2^2) ...
                    / (length(group1)+length(group2)-2));
    d = (mean1 - mean2) / pooledSD;
    stats2(f).cohen_d = d;
end
Parcel2_D = [stats2.cohen_d];
% [minParcel2_D, idxLowest2D] = min(Parcel2_D); 
% display([minParcel2_D, julich_lh.labels.name(idxLowest2D+1)]);

[Parcel2_D_sorted, idxLowest2D_all] = sort(Parcel2_D, 'ascend');
idxLowest2D = idxLowest2D_all(1);  % keep original variable for downstream use
for k = 1:3
    fprintf('Top %d: Cohen''s D = %.4f | Region = %s\n', ...
        k, Parcel2_D_sorted(k), julich_lh.labels.name{idxLowest2D_all(k)+1});
end

% %% 3.1 Plot Parcel with lowest Cohen's D; Economo
% Parcel_brh = zeros(32492,1);
% Parcel_brh(economo_parc(1:end/2, 1) == allTempLbls(idxLowestD)) = -1;
% Parcel_brh(economo_parc(1:end/2, 1) ~= allTempLbls(idxLowestD)) = 0;
% plot_hemispheres(Parcel_brh, inflated_surf_lh, 'views', 'lm'); colormap(flipud(blue));
% set(findall(gcf,'type','axes'),'CLim',[-1 0]);
% 
% %% 3. Plot Parcel with lowest Cohen's D; Julich
% Parcel2_brh = zeros(32492,1);
% Parcel2_brh(julich_lh.cdata(:, 1) == idxLowestD) = -1;
% Parcel2_brh(julich_lh.cdata(:, 1) ~= idxLowestD) = 0;
% plot_hemispheres(Parcel2_brh, inflated_surf_lh, 'views', 'lm'); colormap(flipud(blue));
% set(findall(gcf,'type','axes'),'CLim',[-1 0]);

%% 4. Find the 'Bridgehead'
% Unfolded HippUnfold (used only to get bridgehead vertex indices)
unfold_surf_path = "/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/derivatives/hippunfold_v1.4.1/hippunfold/sub-MNI152NLin2009cAsym/surf/sub-MNI152NLin2009cAsym_hemi-L_space-unfold_den-0p5mm_label-hipp_midthickness.surf.gii";
% HippUnfold surface in MNI space (same vertex count + ordering as unfolded
hipp_mni_surf_path = "/host/verges/tank/data/mo/software/hippunfold-templateflow-main/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_hemi-L_space-T1w_den-0p5mm_label-hipp_inner.surf.gii";
% fsLR-32k surface in MNI space (your “standard neocortical” surface)
fslr_surf_path = "/data_/mica1/01_programs/micapipe-v0.2.0/surfaces/fsLR-32k.L.midthickness.surf.gii";

% Threshold defining "bridgehead" in unfolded space (your example: y < -199.7)
bridgehead_axis = 2;           % MATLAB column index (1=x,2=y,3=z)
bridgehead_thresh = -199.7;    % match your python code
% Contour distances (mm)
medial_dist_mm  =  2.5;
lateral_dist_mm = 15.5;
contour_tol_mm = 0.75; % Contour band tolerance (mm): picks vertices with |d - target| <= tol, % Increase if the band is too sparse/disconnected. 
ap_axis = 2; % Posterior direction definition: In MNI, anterior-posterior is typically Y (2nd coord), anterior = larger Y.

% -------------------------- Load surfaces --------------------------
fprintf("Loading unfolded HippUnfold surface...\n");
[U_V, ~] = load_gifti_surf(unfold_surf_path);
fprintf("Loading MNI-space HippUnfold surface (for mapping to fsLR)...\n");
[H_V, ~] = load_gifti_surf(hipp_mni_surf_path);
fprintf("Loading fsLR-32k MNI-space surface...\n");
[F_V, F_F] = load_gifti_surf(fslr_surf_path);
fprintf("Unfolded mins: [%g %g %g]\n", min(U_V,[],1));
fprintf("Unfolded maxs: [%g %g %g]\n", max(U_V,[],1));

% -------------------- 1) Find bridgehead indices -------------------
bridgehead = find(U_V(:, bridgehead_axis) < bridgehead_thresh);
fprintf("Found %d bridgehead vertices in unfolded space.\n", numel(bridgehead));
if isempty(bridgehead)
    error("No bridgehead vertices found. Check threshold/axis or input surface.");
end

% ---- 2) Map bridgehead vertices to fsLR vertices (MNI space) -------
%   Vertex ordering is consistent between unfolded HippUnfold and MNI-space HippUnfold,
%   so "bridgehead" indices can be used to index H_V.
H_bridge_xyz = H_V(bridgehead, :);

% Map HippUnfold bridgehead xyz -> nearest fsLR vertex
fprintf("Mapping bridgehead points to fsLR vertices via nearest-neighbor...\n");
% knnsearch requires Statistics and Machine Learning Toolbox
fslr_bridge_vidx = knnsearch(F_V, H_bridge_xyz);  % returns fsLR vertex indices (1-based)

% Unique seed set (multi-source distances)
idxVtxSeed = unique(fslr_bridge_vidx);
fprintf("Mapped to %d unique fsLR seed vertices.\n", numel(idxVtxSeed));

bridge_labels = economo_parc(idxVtxSeed);
in_mask = ismember(bridge_labels, allTempLbls);
idxVtxSeed_masked = idxVtxSeed(in_mask);
fprintf("Bridgehead vertices inside TL mask: %d / %d\n", sum(in_mask), numel(idxVtxSeed));

surf_val = zeros(32492,1);
surf_val(idxVtxSeed_masked) = -1;   % <-- changed from idxVtxSeed
obj = plot_hemispheres(surf_val, inflated_surf_lh, 'views', 'lm'); colormap(flipud(blue));
set(findall(gcf,'type','axes'),'CLim',[-1 0]);

%% 5. Calculate the Cohen's D of 'bridgehead'
load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat')

brHeadVal_PTs = mean(Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET(:,idxVtxSeed),2);
brHeadVal_HCs = mean(Voxel_HCs_neofeats_hmz_flip_asym_wPDB.MidPET(:,idxVtxSeed),2);

[~,brHeadP] = ttest2(brHeadVal_PTs, brHeadVal_HCs);
stats(14).ttest_p = brHeadP;
mean1 = mean(brHeadVal_PTs);
mean2 = mean(brHeadVal_HCs);
sd1 = std(brHeadVal_PTs);
sd2 = std(brHeadVal_HCs);
pooledSD = sqrt(((length(brHeadVal_PTs)-1)*sd1^2 + (length(brHeadVal_HCs)-1)*sd2^2) ...
                / (length(brHeadVal_PTs)+length(brHeadVal_HCs)-2));
d = (mean1 - mean2) / pooledSD;
stats(14).cohen_d = d;
display([stats.cohen_d]);

%% 6. Visualization of the group comparison results
% stats(15).cohen_d = stats2(idxLowest2D).cohen_d;
% stats(15).ttest_p = stats2(idxLowest2D).ttest_p;
for k = 1:3
    stats(14+k).cohen_d = stats2(idxLowest2D_all(k)).cohen_d;
    stats(14+k).ttest_p = stats2(idxLowest2D_all(k)).ttest_p;
end

vals = [stats.cohen_d]';
names = ["IJK", "IA", "IB", "HA", "HB", "HC", "TA", "TB", "TC", "TD", "TE", "TF", "TG", "HBI", "EC", "Ia2", "Id9"];
[v_sorted, idx] = sort(vals, 'ascend');
names_sorted = names(idx);

% Map each bar value to a color index (1..256)
vmin = min(v_sorted);
vmax = max(v_sorted);
if vmax == vmin
    cidx = repmat(256, size(v_sorted));
else
    cidx = 256 - round( (v_sorted - vmin) ./ (vmax - vmin) * (256-1) );
end
cidx = max(1, min(256, cidx));

load_paper
ax = subplot(10,7,[24 25 31 32 38 39]);
hold(ax, 'on');
y = 1:numel(v_sorted);
% Draw bars one by one so each has its own color
for i = 1:numel(v_sorted)
    b = barh(y(i), v_sorted(i), 0.72);
    b.FaceColor = blue(cidx(i),:);
    b.EdgeColor = 'none';
end

ax.YDir = 'reverse';
ax.YTick = y;
ax.YTickLabel = names_sorted;
ax.YAxisLocation = 'right';
ax.TickDir = 'out';

% Add value labels at bar ends
xmax = max(v_sorted);
for i = 1:numel(v_sorted)
    text(ax, v_sorted(i) + 0.02*xmax, y(i), sprintf('%.4f', v_sorted(i)), ...
        'VerticalAlignment','middle','HorizontalAlignment','right', ...
        'FontSize',7,'Color',[0.15 0.15 0.15]);
end
grid(ax,'off');
box(ax,'on');
xlim(ax, [-3, 0]);

% Add a colorbar showing Low->High mapping
colormap(ax, blue);
cb = colorbar(ax, 'westoutside');
cb.Ticks = [0 1];
cb.TickLabels = {'Low','High'};
hold(ax,'off');

load_paper
ax = donutRose_SHAP(flipud(v_sorted), flipud(names_sorted), 'Colormap', flipud(blue), 'ShowPercent', false);
% 
% %% 7. save seed index
% save('/host/verges/tank/data/mo/demo/11idxHBI', 'idxVtxSeed');