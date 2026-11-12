%% eFIG6_2: BTTH1 TLE versus BTTH4 frontal-FCD disease controls
% One sequential script; no local functions and no output files are saved.

clear; clc;

%% 1. Load shared data, surfaces, and colors
load_brh
load_demo
load_surface
load_color

DATA_FILE = ...
    '/host/verges/tank/data/mo/demo/04_BTTH1_BTTH4_FrontalFCD_PET.mat';

load(DATA_FILE, ...
    'Parcel_Econo_BTTH1_feats_hmz_flip_asym_wPDB', ...
    'Parcel_Econo_BTTH4_feats_hmz_flip_asym_wPDB', ...
    'Parcel_Julich_BTTH1_feats_hmz_flip_asym_wPDB', ...
    'Parcel_Julich_BTTH4_feats_hmz_flip_asym_wPDB', ...
    'Voxel_BTTH1_neofeats_hmz_flip_asym_wPDB', ...
    'Voxel_BTTH4_neofeats_hmz_flip_asym_wPDB');

%% 2. FIG1E: Economo parcel-wise Cohen's d
% Cohen's d direction: BTTH1 TLE minus BTTH4 frontal FCD.

allTempLbls = [periTempLbls,mesTempLbls,latTempLbls];
numEconoRegions = numel(allTempLbls);
statsEcono = struct();

for f = 1:numEconoRegions
    group1 = Parcel_Econo_BTTH1_feats_hmz_flip_asym_wPDB.MidPET( ...
        :,allTempLbls(f));
    group2 = Parcel_Econo_BTTH4_feats_hmz_flip_asym_wPDB.MidPET( ...
        :,allTempLbls(f));

    [~,p] = ttest2(group1,group2);
    sd1 = std(group1);
    sd2 = std(group2);
    pooledSD = sqrt( ...
        ((numel(group1)-1)*sd1^2+(numel(group2)-1)*sd2^2) / ...
        (numel(group1)+numel(group2)-2));

    statsEcono(f).ttest_p = p;
    statsEcono(f).cohen_d = (mean(group1)-mean(group2))/pooledSD;
end

Econo_D = [statsEcono.cohen_d];
[minEcono_D,idxLowestEcono] = min(Econo_D);
fprintf('\nEconomo minimum Cohen''s d = %.4f; label = %d\n', ...
    minEcono_D,allTempLbls(idxLowestEcono));

%% 3. FIG1E: Julich parcel-wise Cohen's d
numJulichRegions = 148;
statsJulich = struct();

for f = 1:numJulichRegions
    group1 = Parcel_Julich_BTTH1_feats_hmz_flip_asym_wPDB.MidPET(:,f);
    group2 = Parcel_Julich_BTTH4_feats_hmz_flip_asym_wPDB.MidPET(:,f);

    [~,p] = ttest2(group1,group2);
    sd1 = std(group1);
    sd2 = std(group2);
    pooledSD = sqrt( ...
        ((numel(group1)-1)*sd1^2+(numel(group2)-1)*sd2^2) / ...
        (numel(group1)+numel(group2)-2));

    statsJulich(f).ttest_p = p;
    statsJulich(f).cohen_d = (mean(group1)-mean(group2))/pooledSD;
end

Julich_D = [statsJulich.cohen_d];
[Julich_D_sorted,idxLowestJulich_all] = sort(Julich_D,'ascend');

fprintf('\nThree lowest Julich Cohen''s d values:\n');
for k = 1:3
    fprintf('Top %d: d = %.4f | Region = %s\n', ...
        k,Julich_D_sorted(k), ...
        julich_lh.labels.name{idxLowestJulich_all(k)+1});
end

%% 4. FIG1E: Define and display the bridgehead
unfold_surf_path = ...
    "/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/derivatives/hippunfold_v1.4.1/hippunfold/sub-MNI152NLin2009cAsym/surf/sub-MNI152NLin2009cAsym_hemi-L_space-unfold_den-0p5mm_label-hipp_midthickness.surf.gii";
hipp_mni_surf_path = ...
    "/host/verges/tank/data/mo/software/hippunfold-templateflow-main/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_hemi-L_space-T1w_den-0p5mm_label-hipp_inner.surf.gii";
fslr_surf_path = ...
    "/data_/mica1/01_programs/micapipe-v0.2.0/surfaces/fsLR-32k.L.midthickness.surf.gii";

bridgehead_axis = 2;
bridgehead_thresh = -199.7;

[U_V,~] = load_gifti_surf(unfold_surf_path);
[H_V,~] = load_gifti_surf(hipp_mni_surf_path);
[F_V,F_F] = load_gifti_surf(fslr_surf_path);

bridgehead = find(U_V(:,bridgehead_axis)<bridgehead_thresh);
if isempty(bridgehead)
    error('No bridgehead vertices found. Check threshold and axis.');
end

H_bridge_xyz = H_V(bridgehead,:);
fslr_bridge_vidx = knnsearch(F_V,H_bridge_xyz);
idxVtxSeed = unique(fslr_bridge_vidx);

bridge_labels = economo_parc(idxVtxSeed);
in_mask = ismember(bridge_labels,allTempLbls);
idxVtxSeed_masked = idxVtxSeed(in_mask);

fprintf('\nHippUnfold bridgehead vertices: %d\n',numel(bridgehead));
fprintf('Unique fsLR bridgehead vertices: %d\n',numel(idxVtxSeed));
fprintf('Bridgehead vertices inside temporal mask: %d/%d\n', ...
    numel(idxVtxSeed_masked),numel(idxVtxSeed));

surf_val = zeros(32492,1);
surf_val(idxVtxSeed_masked) = -1;
plot_hemispheres(surf_val,inflated_surf_lh,'views','lm');
colormap(flipud(blue));
set(findall(gcf,'type','axes'),'CLim',[-1 0]);

%% 5. FIG1E: Bridgehead Cohen's d
brHeadVal_BTTH1 = mean( ...
    Voxel_BTTH1_neofeats_hmz_flip_asym_wPDB.MidPET(:,idxVtxSeed),2);
brHeadVal_BTTH4 = mean( ...
    Voxel_BTTH4_neofeats_hmz_flip_asym_wPDB.MidPET(:,idxVtxSeed),2);

[~,brHeadP] = ttest2(brHeadVal_BTTH1,brHeadVal_BTTH4);
sd1 = std(brHeadVal_BTTH1);
sd2 = std(brHeadVal_BTTH4);
pooledSD = sqrt( ...
    ((numel(brHeadVal_BTTH1)-1)*sd1^2 + ...
     (numel(brHeadVal_BTTH4)-1)*sd2^2) / ...
    (numel(brHeadVal_BTTH1)+numel(brHeadVal_BTTH4)-2));
brHeadD = (mean(brHeadVal_BTTH1)-mean(brHeadVal_BTTH4))/pooledSD;

fprintf('Bridgehead: d = %.4f; P = %.6g\n',brHeadD,brHeadP);

%% 6. FIG1E: Assemble and display results
econNames = ["IJK","IA","IB","HA","HB","HC", ...
    "TA","TB","TC","TD","TE","TF","TG"];

julichTopNames = strings(1,3);
for k = 1:3
    julichTopNames(k) = string( ...
        julich_lh.labels.name{idxLowestJulich_all(k)+1});
end

names = [econNames,"HBI",julichTopNames];
vals = [Econo_D,brHeadD,Julich_D_sorted(1:3)]';
pvals = [[statsEcono.ttest_p],brHeadP, ...
    [statsJulich(idxLowestJulich_all(1:3)).ttest_p]]';

Results_FIG1E = table(names(:),vals,pvals, ...
    'VariableNames',{'Region', ...
    'CohensD_BTTH1_minus_BTTH4_FrontalFCD','PValue'});
Results_FIG1E = sortrows( ...
    Results_FIG1E,'CohensD_BTTH1_minus_BTTH4_FrontalFCD','ascend');
disp(Results_FIG1E);

[v_sorted,idxSort] = sort(vals,'ascend');
names_sorted = names(idxSort);
vmin = min(v_sorted);
vmax = max(v_sorted);

if vmax==vmin
    cidx = repmat(256,size(v_sorted));
else
    cidx = 256-round((v_sorted-vmin)./(vmax-vmin)*(256-1));
end
cidx = max(1,min(256,cidx));

load_paper
ax = subplot(10,7,[24 25 31 32 38 39]);
hold(ax,'on');
y = 1:numel(v_sorted);

for i = 1:numel(v_sorted)
    b = barh(y(i),v_sorted(i),0.72);
    b.FaceColor = blue(cidx(i),:);
    b.EdgeColor = 'none';
end

ax.YDir = 'reverse';
ax.YTick = y;
ax.YTickLabel = names_sorted;
ax.YAxisLocation = 'right';
ax.TickDir = 'out';

for i = 1:numel(v_sorted)
    text(ax,v_sorted(i),y(i),sprintf(' %.4f',v_sorted(i)), ...
        'VerticalAlignment','middle','HorizontalAlignment','right', ...
        'FontSize',7,'Color',[0.15 0.15 0.15]);
end

grid(ax,'off');
box(ax,'on');
xLower = min(-3,floor(min(v_sorted)*10)/10);
xUpper = max(0,ceil(max(v_sorted)*10)/10);
xlim(ax,[xLower xUpper]);
colormap(ax,blue);
cb = colorbar(ax,'westoutside');
cb.Ticks = [0 1];
cb.TickLabels = {'Low','High'};
hold(ax,'off');

load_paper
donutRose_SHAP(flipud(v_sorted),flipud(names_sorted), ...
    'Colormap',flipud(blue),'ShowPercent',false);

%% 7. FIG1G: Mean cortical PET maps
pet_val_tle = mean( ...
    Voxel_BTTH1_neofeats_hmz_flip_asym_wPDB.MidPET(:,1:32492),1)';
pet_val_fcd = mean( ...
    Voxel_BTTH4_neofeats_hmz_flip_asym_wPDB.MidPET(:,1:32492),1)';

econ_lab = economo_parc(1:end/2,1);
if numel(econ_lab)~=size(F_V,1) || ...
        numel(pet_val_tle)~=size(F_V,1) || ...
        numel(pet_val_fcd)~=size(F_V,1)
    error('Label/PET length mismatch with fsLR vertices.');
end

%% 8. FIG1G: PET values across geodesic contour bands
if isempty(max_dist_mm)
    max_dist_mm = nanmax(dist_roi);
end

levels = 0:contour_step_mm:max_dist_mm;
nL = numel(levels);
mean_tle = nan(nL,1);
ci_tle = nan(nL,1);
mean_fcd = nan(nL,1);
ci_fcd = nan(nL,1);
n_tle = zeros(nL,1);
n_fcd = zeros(nL,1);

for k = 1:nL
    L = levels(k);
    band = is_TL(1:end/2,1) & ...
        abs(dist_roi(1:end/2,1)-L)<=contour_tol_mm;

    band_tle = band & ~isnan(pet_val_tle);
    n_tle(k) = nnz(band_tle);
    if n_tle(k)>1
        bandValues = pet_val_tle(band_tle);
        mean_tle(k) = mean(bandValues,'omitnan');
        ci_tle(k) = 1.96*std(bandValues,'omitnan')/sqrt(n_tle(k));
    end

    band_fcd = band & ~isnan(pet_val_fcd);
    n_fcd(k) = nnz(band_fcd);
    if n_fcd(k)>1
        bandValues = pet_val_fcd(band_fcd);
        mean_fcd(k) = mean(bandValues,'omitnan');
        ci_fcd(k) = 1.96*std(bandValues,'omitnan')/sqrt(n_fcd(k));
    end
end

Results_FIG1G = table( ...
    levels(:),mean_tle,ci_tle,n_tle,mean_fcd,ci_fcd,n_fcd, ...
    'VariableNames',{'distance_mm','mean_PET_BTTH1_TLE', ...
    'CI95_BTTH1_TLE','n_vertices_BTTH1_TLE', ...
    'mean_PET_BTTH4_FrontalFCD','CI95_BTTH4_FrontalFCD', ...
    'n_vertices_BTTH4_FrontalFCD'});
disp(Results_FIG1G(1:min(15,height(Results_FIG1G)),:));

%% 9. FIG1G: Plot mean and 95% CI curves
load_paper
ax_main = subplot(10,7,[44 45 46 47 51 52 53 54 58 59 60 61]);
hold(ax_main,'on');

tle_fill = blue(64,:);
fcd_fill = gray(64,:);
cmap_tle = blue(65:256,:);
cmap_fcd = gray(65:256,:);
nColors = size(cmap_tle,1);

range_tle = max(n_tle)-min(n_tle);
if range_tle==0
    idx_tle = ones(size(n_tle));
else
    idx_tle = round((n_tle-min(n_tle))/range_tle*(nColors-1)+1);
end
idx_tle = max(1,min(nColors,idx_tle));

range_fcd = max(n_fcd)-min(n_fcd);
if range_fcd==0
    idx_fcd = ones(size(n_fcd));
else
    idx_fcd = round((n_fcd-min(n_fcd))/range_fcd*(nColors-1)+1);
end
idx_fcd = max(1,min(nColors,idx_fcd));

valid_tle = ~isnan(mean_tle) & ~isnan(ci_tle);
valid_fcd = ~isnan(mean_fcd) & ~isnan(ci_fcd);
x_tle = levels(valid_tle)';
x_fcd = levels(valid_fcd)';

fill(ax_main,[x_fcd;flipud(x_fcd)], ...
    [mean_fcd(valid_fcd)+ci_fcd(valid_fcd); ...
    flipud(mean_fcd(valid_fcd)-ci_fcd(valid_fcd))], ...
    fcd_fill,'EdgeColor','none','FaceAlpha',0.3);

fill(ax_main,[x_tle;flipud(x_tle)], ...
    [mean_tle(valid_tle)+ci_tle(valid_tle); ...
    flipud(mean_tle(valid_tle)-ci_tle(valid_tle))], ...
    tle_fill,'EdgeColor','none','FaceAlpha',0.3);

lev_fcd = levels(valid_fcd);
m_fcd = mean_fcd(valid_fcd);
idx_fcd_valid = idx_fcd(valid_fcd);
for i = 1:numel(lev_fcd)-1
    plot(ax_main,lev_fcd(i:i+1),m_fcd(i:i+1),'-', ...
        'Color',cmap_fcd(idx_fcd_valid(i),:),'LineWidth',2);
end

lev_tle = levels(valid_tle);
m_tle = mean_tle(valid_tle);
idx_tle_valid = idx_tle(valid_tle);
for i = 1:numel(lev_tle)-1
    plot(ax_main,lev_tle(i:i+1),m_tle(i:i+1),'-', ...
        'Color',cmap_tle(idx_tle_valid(i),:),'LineWidth',2);
end

box(ax_main,'off');
xlabel(ax_main,'Geodesic distance to bridgehead');
ylabel(ax_main,'Mean PET within contour band');
ylim(ax_main,[-4.5 0.5]);
xlim(ax_main,[-5 130]);

pos = ax_main.Position;
ax_cb_tle = axes('Position',[pos(1)+pos(3)+0.01,pos(2),0.015,pos(4)]);
colormap(ax_cb_tle,cmap_tle);
cb_tle = colorbar(ax_cb_tle,'Location','eastoutside');
ax_cb_tle.Visible = 'off';
cb_tle.Label.String = 'n_{vertices,TLE}';
cb_tle.Ticks = [0 1];
cb_tle.TickLabels = {num2str(min(n_tle)),num2str(max(n_tle))};

ax_cb_fcd = axes('Position',[pos(1)+pos(3)+0.01,pos(2),0.015,pos(4)]);
colormap(ax_cb_fcd,cmap_fcd);
cb_fcd = colorbar(ax_cb_fcd,'Location','eastoutside');
ax_cb_fcd.Visible = 'off';
cb_fcd.Position(1) = cb_tle.Position(1)+cb_tle.Position(3)+0.02;
cb_fcd.Label.String = 'n_{vertices,frontal FCD}';
cb_fcd.Ticks = [0 1];
cb_fcd.TickLabels = {num2str(min(n_fcd)),num2str(max(n_fcd))};

%% 10. FIG1G: Distance-gradient permutation tests
nSpin = min(123,numel(levels));

fprintf('\nBTTH4 frontal-FCD distance-gradient permutation test:\n');
spin_permutation_test( ...
    levels(1:nSpin),mean_fcd(1:nSpin),'n_perm',10000);

fprintf('\nBTTH1 TLE distance-gradient permutation test:\n');
spin_permutation_test( ...
    levels(1:nSpin),mean_tle(1:nSpin),'n_perm',10000);

fprintf('\nCompleted frontal-FCD disease-control analyses.\n');
