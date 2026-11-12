%% eFIG6_12_1: TLE ipsilateral PET asymmetry referenced to non-temporal non-insular FCD right-minus-left PET asymmetry
% Input from Step0_12: BTTH1 TLE = ipsilateral - contralateral asymmetry w-scores;
% BTTH4 non-temporal non-insular FCD reference = right - left asymmetry (no lesion-side flipping).
% Contains: FIG1E/FIG1G (no surface plot, no rose plot), FIG2A cyto layers,
% FIG3A MidPET (midsurface PET SUVR) comparison, FIG4A HippoPET main effect,
% FIG7B surgical metrics / logistic regression using the new w-scores.
% One sequential script; no local functions and no output files are saved.

clear; clc; rng('default');

%% 1. Load shared data, surfaces, and colors
load_brh
load_demo
load_surface
load_color

nVhemi = 32492;

DATA_FILE = ...
    '/host/verges/tank/data/mo/demo/12_BTTH1_BTTH4_nonTemporalNonInsularFCD_PET_RL.mat';

load(DATA_FILE, ...
    'Parcel_Econo_BTTH1_feats_hmz_flip_asym_wPDB', ...
    'Parcel_Julich_BTTH1_feats_hmz_flip_asym_wPDB', ...
    'Voxel_BTTH1_neofeats_hmz_flip_asym_wPDB', ...
    'Voxel_BTTH4_neofeats_hmz_flip_asym_wPDB');

% Align the BTTH1 rows of the new file with Demo_disPTs rows.
% Step0_12 keeps only patients passing inclusion, in Demo_disPTs order.
% If ID_BTTH1 was also saved in the .mat file it is used; otherwise the
% number of rows must equal the number of patients in Demo_disPTs.
fileVars = who('-file',DATA_FILE);
nBTTH1 = size(Voxel_BTTH1_neofeats_hmz_flip_asym_wPDB.MidPET,1);
nPT_all = numel(Demo_disPTs.ID);
if ismember('ID_BTTH1',fileVars)
    S_id = load(DATA_FILE,'ID_BTTH1');
    ptIdx = find(ismember(string(Demo_disPTs.ID(:)),string(S_id.ID_BTTH1(:))));
elseif nBTTH1==nPT_all
    ptIdx = (1:nPT_all)';
else
    error(['BTTH1 rows (%d) differ from Demo_disPTs rows (%d) and ID_BTTH1 ' ...
        'is not in the .mat file. Add ID_BTTH1 to the save list of Step0_12.'], ...
        nBTTH1,nPT_all);
end
if numel(ptIdx)~=nBTTH1
    error('Could not align BTTH1 rows with Demo_disPTs (%d vs %d).', ...
        numel(ptIdx),nBTTH1);
end
fprintf('\nBTTH1 TLE patients: %d; BTTH4 reference patients: %d\n', ...
    nBTTH1,size(Voxel_BTTH4_neofeats_hmz_flip_asym_wPDB.MidPET,1));

%% 2. FIG1E: Economo parcel-wise one-sample Cohen's d versus zero
allTempLbls = [periTempLbls,mesTempLbls,latTempLbls];
numEconoRegions = numel(allTempLbls);
statsEcono = struct();

for f = 1:numEconoRegions
    group1 = Parcel_Econo_BTTH1_feats_hmz_flip_asym_wPDB.MidPET( ...
        :,allTempLbls(f));
    [~,p] = ttest(group1,0);

    statsEcono(f).ttest_p = p;
    statsEcono(f).cohen_d = mean(group1,'omitnan')/std(group1,'omitnan');
end

Econo_D = [statsEcono.cohen_d];
[minEcono_D,idxLowestEcono] = min(Econo_D);
fprintf('\nEconomo minimum Cohen''s d = %.4f; label = %d\n', ...
    minEcono_D,allTempLbls(idxLowestEcono));

%% 3. FIG1E: Julich parcel-wise one-sample Cohen's d versus zero
numJulichRegions = 148;
statsJulich = struct();

for f = 1:numJulichRegions
    group1 = Parcel_Julich_BTTH1_feats_hmz_flip_asym_wPDB.MidPET(:,f);
    [~,p] = ttest(group1,0);

    statsJulich(f).ttest_p = p;
    statsJulich(f).cohen_d = mean(group1,'omitnan')/std(group1,'omitnan');
end

Julich_D = [statsJulich.cohen_d];
[Julich_D_sorted,idxLowestJulich_all] = sort(Julich_D,'ascend');

fprintf('\nThree lowest Julich Cohen''s d values:\n');
for k = 1:3
    fprintf('Top %d: d = %.4f | Region = %s\n', ...
        k,Julich_D_sorted(k), ...
        julich_lh.labels.name{idxLowestJulich_all(k)+1});
end

%% 4. FIG1E: Define the bridgehead (no surface display)
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

%% 5. FIG1E: Bridgehead one-sample Cohen's d versus zero
brHeadVal_BTTH1 = mean( ...
    Voxel_BTTH1_neofeats_hmz_flip_asym_wPDB.MidPET(:,idxVtxSeed),2);
[~,brHeadP] = ttest(brHeadVal_BTTH1,0);
brHeadD = mean(brHeadVal_BTTH1,'omitnan')/std(brHeadVal_BTTH1,'omitnan');

fprintf('Bridgehead: d = %.4f; P = %.6g\n',brHeadD,brHeadP);

%% 6. FIG1E: Assemble and display results (bar plot only)
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
    'VariableNames',{'Region','CohensD_BTTH1_vs_zero','PValue'});
Results_FIG1E = sortrows( ...
    Results_FIG1E,'CohensD_BTTH1_vs_zero','ascend');
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
figure('Units','centimeters','Position',[2 2 21 29.7]);
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

%% 7. FIG1G: Mean TLE cortical PET map
pet_val_tle = mean( ...
    Voxel_BTTH1_neofeats_hmz_flip_asym_wPDB.MidPET(:,1:32492),1)';

econ_lab = economo_parc(1:end/2,1);
if numel(econ_lab)~=size(F_V,1) || ...
        numel(pet_val_tle)~=size(F_V,1)
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
n_tle = zeros(nL,1);
% Because the non-temporal non-insular FCD right-minus-left asymmetry defines the
% normative model, its expected w-score is zero at every distance.
mean_non_temporal_fcd_control = zeros(nL,1);

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
end

Results_FIG1G = table( ...
    levels(:),mean_tle,ci_tle,n_tle,mean_non_temporal_fcd_control, ...
    'VariableNames',{'distance_mm','mean_PET_BTTH1_TLE', ...
    'CI95_BTTH1_TLE','n_vertices_BTTH1_TLE', ...
    'mean_PET_nonTemporalNonInsularFCD_RL_asym_reference'});
disp(Results_FIG1G(1:min(15,height(Results_FIG1G)),:));

%% 9. FIG1G: Plot mean and 95% CI curves (same figure as the bar plot)
ax_main = subplot(10,7,[44 45 46 47 51 52 53 54 58 59 60 61]);
hold(ax_main,'on');

tle_fill = blue(64,:);
control_color = [0.55 0.55 0.55];
cmap_tle = blue(65:256,:);
nColors = size(cmap_tle,1);

range_tle = max(n_tle)-min(n_tle);
if range_tle==0
    idx_tle = ones(size(n_tle));
else
    idx_tle = round((n_tle-min(n_tle))/range_tle*(nColors-1)+1);
end
idx_tle = max(1,min(nColors,idx_tle));

valid_tle = ~isnan(mean_tle) & ~isnan(ci_tle);
x_tle = levels(valid_tle)';

h_control = plot(ax_main,levels,mean_non_temporal_fcd_control,'-', ...
    'Color',control_color,'LineWidth',2);

fill(ax_main,[x_tle;flipud(x_tle)], ...
    [mean_tle(valid_tle)+ci_tle(valid_tle); ...
    flipud(mean_tle(valid_tle)-ci_tle(valid_tle))], ...
    tle_fill,'EdgeColor','none','FaceAlpha',0.3);

lev_tle = levels(valid_tle);
m_tle = mean_tle(valid_tle);
idx_tle_valid = idx_tle(valid_tle);
for i = 1:numel(lev_tle)-1
    plot(ax_main,lev_tle(i:i+1),m_tle(i:i+1),'-', ...
        'Color',cmap_tle(idx_tle_valid(i),:),'LineWidth',2);
end

h_tle = plot(ax_main,nan,nan,'-','Color',blue(192,:),'LineWidth',2);
legend(ax_main,[h_tle,h_control], ...
    {'TLE ipsilateral-contralateral asymmetry', ...
    'Non-temporal non-insular FCD R-L asymmetry reference'}, ...
    'Location','southeast','Box','off');

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

%% 10. FIG1G: Distance-gradient permutation tests
nSpin = min(123,numel(levels));

fprintf(['\nBTTH1 TLE ipsilateral-contralateral PET asymmetry, referenced to ' ...
    'non-temporal non-insular FCD right-minus-left PET asymmetry:\n']);
spin_permutation_test( ...
    levels(1:nSpin),mean_tle(1:nSpin),'n_perm',10000);

%% 11. FIG2A: Economo cytoarchitecture-layer comparison
% Subject-wise mean w-score within each of the 5 layer classes (temporal
% mask, left/ipsilateral half). Two-sample test: FCD reference vs TLE.
condition_names = {'1','2','3','4','5'};
num_cond = numel(condition_names);

econ_L   = economo_parc(1:nVhemi);
layers_L = economo_parc_5layers(1:nVhemi);
W_TLE_L  = Voxel_BTTH1_neofeats_hmz_flip_asym_wPDB.MidPET(:,1:nVhemi);
W_REF_L  = Voxel_BTTH4_neofeats_hmz_flip_asym_wPDB.MidPET(:,1:nVhemi);

cyto_REF = zeros(size(W_REF_L,1),num_cond);
cyto_TLE = zeros(size(W_TLE_L,1),num_cond);
nVtx_layer = zeros(1,num_cond);
for i = 1:num_cond
    idxTmp = ismember(econ_L,allTempLbls) & (layers_L == i);
    nVtx_layer(i) = nnz(idxTmp);
    cyto_REF(:,i) = mean(W_REF_L(:,idxTmp),2,'omitnan');
    cyto_TLE(:,i) = mean(W_TLE_L(:,idxTmp),2,'omitnan');
end

cyto_d = zeros(1,num_cond);
cyto_t = zeros(1,num_cond);
cyto_p = zeros(1,num_cond);
for i = 1:num_cond
    hc = cyto_REF(:,i); hc = hc(~isnan(hc));
    pt = cyto_TLE(:,i); pt = pt(~isnan(pt));
    [~,p,~,st] = ttest2(hc,pt);
    cyto_t(i) = st.tstat;
    cyto_p(i) = p;
    n_hc = numel(hc); n_pt = numel(pt);
    pooled_std = sqrt(((n_hc-1)*std(hc)^2+(n_pt-1)*std(pt)^2)/(n_hc+n_pt-2));
    cyto_d(i) = (mean(pt)-mean(hc))/pooled_std;
end
[~,~,~,cyto_adj_p] = fdr_bh(cyto_p,0.05,'pdep','yes');

Results_Cyto = table(string(condition_names(:)),nVtx_layer(:), ...
    cyto_d(:),cyto_t(:),cyto_p(:),cyto_adj_p(:), ...
    'VariableNames',{'Layer','nVertices','CohensD_TLE_minus_ref', ...
    'tstat','PValue','PValue_FDR'});
fprintf('\nFIG2A cytoarchitecture-layer comparison (FCD reference vs TLE):\n');
disp(Results_Cyto);

violin_offset = 0.22; box_offset = 0.12; violin_width_scale = 0.20;
ref_alpha = 0.20; tle_alpha = 1.00; lw = 0.1; box_half_w = 0.04;
cond_x = 0:num_cond-1;
load_paper
figure('Units','centimeters','Position',[2 2 21 29.7]);
subplot(10,7,[1 2 3 8 9 10]); hold on;
for c = 1:num_cond
    for g = 1:2
        if g==1
            dd = cyto_REF(:,c); sgn = -1; fa = ref_alpha;
        else
            dd = cyto_TLE(:,c); sgn = 1;  fa = tle_alpha;
        end
        dd = dd(~isnan(dd));
        col = colors_cyto_nor(c,:);
        vx = cond_x(c)+sgn*violin_offset;
        [fk,xi] = ksdensity(dd,'Bandwidth',0.3);
        fk = violin_width_scale*(fk/max(fk));
        fill([vx+sgn*fk, vx*ones(1,numel(xi))],[xi, fliplr(xi)],col, ...
            'FaceAlpha',fa,'EdgeColor','k','LineWidth',lw);
        bx = cond_x(c)+sgn*box_offset;
        q = prctile(dd,[25 50 75]); iq = q(3)-q(1);
        wlo = min(dd(dd >= q(1)-0.5*iq)); whi = max(dd(dd <= q(3)+0.5*iq));
        fill([bx-box_half_w,bx+box_half_w,bx+box_half_w,bx-box_half_w], ...
            [q(1),q(1),q(3),q(3)],col,'FaceAlpha',fa,'EdgeColor','k','LineWidth',lw);
        plot([bx-box_half_w,bx+box_half_w],[q(2),q(2)],'-','Color','w','LineWidth',lw);
        plot([bx bx],[wlo q(1)],'k-','LineWidth',lw);
        plot([bx bx],[q(3) whi],'k-','LineWidth',lw);
    end
end
med_ref = median(cyto_REF,'omitnan'); med_tle = median(cyto_TLE,'omitnan');
for c = 1:num_cond-1
    col_avg = (colors_cyto_nor(c,:)+colors_cyto_nor(c+1,:))/2;
    plot(cond_x([c c+1])-box_offset+box_half_w*0.5,med_ref([c c+1]),'--','Color',col_avg,'LineWidth',lw);
    plot(cond_x([c c+1])+box_offset-box_half_w*0.5,med_tle([c c+1]),'--','Color',col_avg,'LineWidth',lw);
end
set(gca,'XTick',cond_x,'XTickLabel',condition_names,'FontSize',5,'Box','off');
xlim([cond_x(1)-violin_offset-violin_width_scale-0.08, ...
      cond_x(end)+violin_offset+violin_width_scale+0.08]);
ylabel('Mean w-score (left: FCD reference, right: TLE)');

%% 12. FIG3A: Midsurface PET SUVR (MidPET) within-group comparison
% Vertex-wise group-mean maps inside the temporal mask. Left = TLE (BTTH1),
% right = FCD reference (BTTH4). Two-sample t-test and pooled Cohen's d.
inMaskL = ismember(economo_parc(1:nVhemi),allTempLbls);
d_h = mean(Voxel_BTTH1_neofeats_hmz_flip_asym_wPDB.MidPET(:,1:nVhemi),1)';
d_d = mean(Voxel_BTTH4_neofeats_hmz_flip_asym_wPDB.MidPET(:,1:nVhemi),1)';
d_h = d_h(inMaskL); d_d = d_d(inMaskL);
d_h = d_h(~isoutlier(d_h,'gesd'));
d_d = d_d(~isoutlier(d_d,'gesd'));

[~,p3] = ttest2(d_h,d_d);
n1 = numel(d_h); n2 = numel(d_d);
pooledSD = sqrt(((n1-1)*std(d_h)^2+(n2-1)*std(d_d)^2)/(n1+n2-2));
d3 = (mean(d_h)-mean(d_d))/pooledSD;
fprintf('\nFIG3A MidPET: TLE vs FCD reference: p = %.3e, d = %.3f (n vertices = %d vs %d)\n', ...
    p3,d3,n1,n2);

left_col = blue(200,:); right_col = gray(200,:);
cx = 1; v_off = 0.2; s_off = 0.2; box_w = 0.08; v_w = 0.15; lw3 = 0.5;
load_paper
figure('Units','centimeters','Position',[2 2 21 29.7]);
subplot(10,7,35); hold on;
for g = 1:2
    if g==1
        dd = d_h; col = left_col; sgn = -1; fa = 1.0;
    else
        dd = d_d; col = right_col; sgn = 1; fa = 0.30;
    end
    vx = cx+sgn*v_off;
    [fk,xi] = ksdensity(dd,'Bandwidth',0.3); fk = v_w*(fk/max(fk));
    fill([vx+sgn*fk, vx*ones(1,numel(xi))],[xi, fliplr(xi)],col, ...
        'FaceAlpha',fa,'EdgeColor','k','LineWidth',lw3,'EdgeAlpha',1);
    bx = cx+sgn*s_off;
    q = prctile(dd,[25 50 75]); iq = q(3)-q(1);
    wlo = min(dd(dd >= q(1)-0.5*iq)); whi = max(dd(dd <= q(3)+0.5*iq));
    fill([bx-box_w bx+box_w bx+box_w bx-box_w],[q(1) q(1) q(3) q(3)],col, ...
        'FaceAlpha',fa,'EdgeColor','k','LineWidth',lw3,'EdgeAlpha',1);
    plot([bx-box_w bx+box_w],[q(2) q(2)],'-','Color','w','LineWidth',lw3);
    plot([bx bx],[wlo q(1)],'k-','LineWidth',lw3);
    plot([bx bx],[q(3) whi],'k-','LineWidth',lw3);
end
set(gca,'XTick',[cx-v_off cx+v_off],'XTickLabel',{'TLE','FCD ref'}, ...
    'FontSize',5,'Box','off');
xlim([cx-v_off-v_w-0.08, cx+v_off+v_w+0.08]);
text(0.02,0.95,sprintf('p=%.3e, d=%.3f',p3,d3),'Units','normalized', ...
    'FontSize',6,'VerticalAlignment','top');
title('MidPET (midsurface PET SUVR)','FontSize',6);

%% 13. FIG4A: HippoPET main effect on neocortical MidPET (no surface plots)
% HippoPET is not recomputed in Step0_12; it is taken from the 02 file
% (hippocampal w-score asymmetry), aligned to the patients in this analysis.
S02 = load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat', ...
    'Voxel_PTs_subfeats_hmz_flip_asym_wPDB');
HippoPET_PTs = double(S02.Voxel_PTs_subfeats_hmz_flip_asym_wPDB.HippoPET(ptIdx,:));
hippoPET_mean = mean(HippoPET_PTs,2);

Age_PT = double(Demo_disPTs.AgeSurg(ptIdx));
Sex_PT = double(Sex_all2_num(length(Demo_HCs.ID)+ptIdx));
Dur_PT = double(Demo_disPTs.Duration(ptIdx));
SZFreq_raw = Demo_disPTs.SZFreq(ptIdx);
SZFreq_PT = zeros(size(SZFreq_raw));
SZFreq_PT(strcmp(SZFreq_raw,'Daily'))   = 1;
SZFreq_PT(strcmp(SZFreq_raw,'Weekly'))  = 2;
SZFreq_PT(strcmp(SZFreq_raw,'Monthly')) = 3;
SZFreq_PT(strcmp(SZFreq_raw,'Yearly'))  = 4;

Age_Cov = term(Age_PT); Sex_Cov = term(Sex_PT);
Dur_Cov = term(Dur_PT); SZFreq_Cov = term(SZFreq_PT);

PET_lh = Voxel_BTTH1_neofeats_hmz_flip_asym_wPDB.MidPET(:,1:nVhemi);
medialIdx = find(economo_parc(1:nVhemi) ~= 0);
PET_lh(:,setdiff(1:nVhemi,medialIdx)) = 0;   % zero the medial wall

Hippo_Var = term(hippoPET_mean);
Model = 1 + Hippo_Var + Age_Cov + Sex_Cov + Dur_Cov + SZFreq_Cov;
slm_hippo = SurfStatLinMod(PET_lh,Model);
slm_hippo_Pos = SurfStatT(slm_hippo, hippoPET_mean);
slm_hippo_Neg = SurfStatT(slm_hippo,-hippoPET_mean);

p_Pos = 1-tcdf(slm_hippo_Pos.t(:),slm_hippo_Pos.df);
p_Neg = 1-tcdf(slm_hippo_Neg.t(:),slm_hippo_Neg.df);
validV = (economo_parc(1:nVhemi)~=0) & isfinite(p_Pos) & isfinite(p_Neg);

h_Pos = false(nVhemi,1); h_Neg = false(nVhemi,1);
[hp,~] = fdr_bh(p_Pos(validV),0.05); h_Pos(validV) = hp;
[hn,~] = fdr_bh(p_Neg(validV),0.05); h_Neg(validV) = hn;

D_Neg = slm_hippo_Neg.ef(:);
D_Pos = slm_hippo_Pos.ef(:);
inTemp = ismember(economo_parc(1:nVhemi),allTempLbls) & validV;

tail = ["Negative (lower hippoPET -> lower neocortical PET w)"; ...
        "Positive"];
nSigUncorr_all  = [nnz(p_Neg<0.05 & validV); nnz(p_Pos<0.05 & validV)];
nSigUncorr_temp = [nnz(p_Neg<0.05 & inTemp); nnz(p_Pos<0.05 & inTemp)];
nSigFDR_all     = [nnz(h_Neg); nnz(h_Pos)];
nSigFDR_temp    = [nnz(h_Neg & inTemp); nnz(h_Pos & inTemp)];
meanD_temp      = [mean(D_Neg(inTemp),'omitnan'); mean(D_Pos(inTemp),'omitnan')];
maxD_temp       = [max(D_Neg(inTemp)); max(D_Pos(inTemp))];

Results_FIG4A = table(tail,nSigUncorr_all,nSigUncorr_temp, ...
    nSigFDR_all,nSigFDR_temp,meanD_temp,maxD_temp, ...
    'VariableNames',{'Tail','nSig_p05_all','nSig_p05_temporal', ...
    'nSig_FDR_all','nSig_FDR_temporal','MeanEffect_temporal','MaxEffect_temporal'});
fprintf('\nFIG4A HippoPET main effect on neocortical MidPET (n = %d):\n',nBTTH1);
disp(Results_FIG4A);

[r4,pr4] = corr(hippoPET_mean,mean(PET_lh(:,inTemp),2), ...
    'Rows','complete');
fprintf('Pearson r (mean HippoPET vs mean temporal MidPET w-score) = %.3f, p = %.3e\n', ...
    r4,pr4);

% HippoPET distribution (left hippocampal half, GESD outliers removed)
featData = double(mean(HippoPET_PTs,1)');
featData = featData(1:end/2,1);
featData = featData(~isoutlier(featData,'gesd'));
n_scatter = 500; bw = 0.12; x_pos = 1;
color_g = cmap_yellow(180,:); solid_c = [0 0 0];
load_paper
figure('Units','centimeters','Position',[2 2 21 29.7]);
ax = subplot(10,7,35); hold(ax,'on');
q1 = prctile(featData,25); q3 = prctile(featData,75); mn = mean(featData);
iqr_v = q3-q1;
w_lo = max(featData(featData >= q1-1.5*iqr_v));
w_hi = min(featData(featData <= q3+1.5*iqr_v));
n_plot = min(n_scatter,numel(featData));
s_idx = randperm(numel(featData),n_plot);
jitter = (rand(n_plot,1)-0.5)*bw*1.8;
scatter(ax,x_pos+jitter,featData(s_idx),4,color_g,'filled', ...
    'MarkerFaceAlpha',0.25,'MarkerEdgeColor','none');
fill(ax,[x_pos-bw,x_pos+bw,x_pos+bw,x_pos-bw],[q1,q1,q3,q3],color_g, ...
    'EdgeColor',solid_c,'LineWidth',0.5,'FaceAlpha',0.75);
plot(ax,[x_pos-bw,x_pos+bw],[mn,mn],'-','Color',solid_c,'LineWidth',2);
plot(ax,[x_pos,x_pos],[w_lo,q1],'-','Color',solid_c,'LineWidth',0.5);
plot(ax,[x_pos,x_pos],[q3,w_hi],'-','Color',solid_c,'LineWidth',0.5);
plot(ax,[x_pos-bw*0.5,x_pos+bw*0.5],[w_lo,w_lo],'-','Color',solid_c,'LineWidth',0.5);
plot(ax,[x_pos-bw*0.5,x_pos+bw*0.5],[w_hi,w_hi],'-','Color',solid_c,'LineWidth',0.5);
pad = max((max(featData)-min(featData))*0.15,0.2);
box(ax,'off'); xlim(ax,[0.6 1.4]);
ylim(ax,[min(featData)-pad, max(featData)+pad]);
xticks(ax,[]); ylabel(ax,'HippoPET','Interpreter','none','FontSize',8);
set(ax,'FontSize',7);

%% 14. FIG7B: Surgical metrics with the new MidPET w-scores
N = nPT_all;
surgVol = nan(N,1);
surgArea = nan(N,1);

petMean = nan(N,1); petLCI = nan(N,1); petUCI = nan(N,1);
petP025 = nan(N,1); petP975 = nan(N,1); petN = nan(N,1);

wmapMean = nan(N,1); wmapLCI = nan(N,1); wmapUCI = nan(N,1);
wmapP025 = nan(N,1); wmapP975 = nan(N,1);

wmapN = nan(N,1);
wmapVolSum = nan(N,1);
wmapVolSurgSum = nan(N,1);
wmapVolSurgProp = nan(N,1);
wmapVolSurgN = nan(N,1);

petP975_ExtAll = nan(N,1);
petP975_ExtTemp = nan(N,1);
wmapP975_ExtAll = nan(N,1);
wmapP975_ExtTemp = nan(N,1);
ratio_P975_ExtAll = nan(N,1);
ratio_P975_ExtTemp = nan(N,1);

maskAll = double(ismember(economo_parc,allTempLbls));
maskTemp = double(ismember(economo_parc,pureTempLbls));
PATH_PTs = '/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/derivatives/';

for k = 1:nBTTH1
    i = ptIdx(k);   % row in Demo_disPTs

    % individual w-score PET (new reference)
    indi_MidPET = Voxel_BTTH1_neofeats_hmz_flip_asym_wPDB.MidPET(k,1:32492)';

    % abnormal w-score map
    indi_MidPET_ThrVal = indi_MidPET;
    indi_MidPET_ThrVal(~(maskAll(1:end/2,1))) = 0;
    indi_MidPET_ThrVal(indi_MidPET_ThrVal > -1.96) = 0;

    % surgical cavity
    File_LH = strcat(PATH_PTs,'surgMask/',Demo_disPTs.ID{i},'/ses-post/',Demo_disPTs.ID{i},'_ses-post_hemi-L_space-nativepro_surf-fsLR-32k_label-cavity.func.gii');
    File_RH = strcat(PATH_PTs,'surgMask/',Demo_disPTs.ID{i},'/ses-post/',Demo_disPTs.ID{i},'_ses-post_hemi-R_space-nativepro_surf-fsLR-32k_label-cavity.func.gii');
    assume_surg_raw = [gifti(File_LH).cdata; gifti(File_RH).cdata];

    if strcmp(Demo_disPTs.Side{i},'R')
        assume_surg_temp = [assume_surg_raw(end/2+1:end,1); assume_surg_raw(1:end/2,1)];
    else
        assume_surg_temp = assume_surg_raw;
    end
    in_mask = ismember(economo_parc,allTempLbls);
    assume_surg_temp(~in_mask) = 0;
    assume_surg = assume_surg_temp(1:end/2,1);

    % 1) Resection size within TL
    idxTempSurg = (assume_surg) & maskAll(1:end/2,1);
    surgVol(i) = sum(idxTempSurg);

    % 2) PET within resection (TL) + 95% CI
    petVals = indi_MidPET(idxTempSurg,1);
    petVals = petVals(~isnan(petVals));
    petN(i) = numel(petVals);
    if petN(i) >= 2
        mu = mean(petVals); sd = std(petVals,0);
        half = tinv(0.975,petN(i)-1)*sd/sqrt(petN(i));
        petMean(i) = mu; petLCI(i) = mu-half; petUCI(i) = mu+half;
        pr = prctile(petVals,[2.5 97.5]);
        petP025(i) = pr(1); petP975(i) = pr(2);
    elseif petN(i) == 1
        petMean(i) = petVals; petP025(i) = petVals; petP975(i) = petVals;
    end

    % 3) Proportion of abnormal w-score resected (TL)
    wmapThr = indi_MidPET_ThrVal(:,1) <= -1.96;
    wmapN(i) = sum(wmapThr);
    idxWmapTempSurg = wmapThr & idxTempSurg;
    wmapVolSum(i) = sum(wmapThr);
    wmapVolSurgSum(i) = sum(idxWmapTempSurg);
    if wmapVolSum(i) > 0
        wmapVolSurgProp(i) = wmapVolSurgSum(i)/wmapVolSum(i);
    else
        wmapVolSurgProp(i) = 0;
    end

    % 4) w-score within resection (TL) + 95% CI
    wmapVolSurg = indi_MidPET_ThrVal(idxWmapTempSurg,1);
    wmapVolSurg = wmapVolSurg(~isnan(wmapVolSurg));
    wmapVolSurgN(i) = numel(wmapVolSurg);
    if wmapVolSurgN(i) >= 2
        mu = mean(wmapVolSurg); sd = std(wmapVolSurg,0);
        half = tinv(0.975,wmapVolSurgN(i)-1)*sd/sqrt(wmapVolSurgN(i));
        wmapMean(i) = mu; wmapLCI(i) = mu-half; wmapUCI(i) = mu+half;
        pr = prctile(wmapVolSurg,[2.5 97.5]);
        wmapP025(i) = pr(1); wmapP975(i) = pr(2);
    elseif wmapVolSurgN(i) == 1
        wmapMean(i) = wmapVolSurg; wmapP025(i) = wmapVolSurg; wmapP975(i) = wmapVolSurg;
    end

    petExtAllMin = min(indi_MidPET(maskAll(1:end/2,1) == 1,1));
    petExtAllMax = max(indi_MidPET(maskAll(1:end/2,1) == 1,1));
    petExtTempMin = min(indi_MidPET(maskTemp(1:end/2,1) == 1,1));
    petExtTempMax = max(indi_MidPET(maskTemp(1:end/2,1) == 1,1));

    petP975_ExtAll(i)  = (petP975(i)-petExtAllMin)/(petExtAllMax-petExtAllMin);
    petP975_ExtTemp(i) = (petP975(i)-petExtTempMin)/(petExtTempMax-petExtTempMin);
    ratio_P975_ExtAll(i)  = petP975(i)/petExtAllMax;
    ratio_P975_ExtTemp(i) = petP975(i)/petExtTempMax;
end
surgVolProp = surgVol/sum(maskAll(1:end/2,1));

% Hippocampal volume asymmetry (independent of PET reference)
hipp_data = readtable('/host/verges/tank/data/mo/demo/52all_hippocampal_volumes_normalized.csv');
normHippoVol_IPS  = nan(N,1);
normHippoVol_CON  = nan(N,1);
normHippoVol_asym = nan(N,1);
for i = 1:N
    subj = Demo_disPTs.ID{i};
    side = Demo_disPTs.Side{i};
    if strcmp(side,'L'), contra_side = 'R'; else, contra_side = 'L'; end
    subj_rows = strcmp(hipp_data.subject,subj);
    if sum(subj_rows) == 0
        warning('Subject %s not found in hippocampal volume file',subj);
        continue
    end
    ipsi_row   = subj_rows & strcmp(hipp_data.hemi,side);
    contra_row = subj_rows & strcmp(hipp_data.hemi,contra_side);
    if sum(ipsi_row) == 0 || sum(contra_row) == 0
        warning('Missing hemisphere for subject %s',subj);
        continue
    end
    normHippoVol_IPS(i) = hipp_data.normalized_hipp_vol(ipsi_row);
    normHippoVol_CON(i) = hipp_data.normalized_hipp_vol(contra_row);
    normHippoVol_asym(i) = (normHippoVol_IPS(i)-normHippoVol_CON(i)) / ...
        (0.5*(abs(normHippoVol_IPS(i))+abs(normHippoVol_CON(i))));
end

T = table((1:N)',surgVol,surgArea,surgVolProp, ...
    petMean,petLCI,petUCI,petP025,petP975, ...
    wmapN,wmapVolSurgN,wmapVolSurgProp, ...
    wmapMean,wmapLCI,wmapUCI,wmapP025,wmapP975, ...
    petP975_ExtAll,petP975_ExtTemp,wmapP975_ExtAll,wmapP975_ExtTemp,normHippoVol_asym, ...
    'VariableNames',{'Subj','surgVol','surgArea','surgVolProp', ...
    'petMean','petLCI','petUCI','PET_p2p5','petP975', ...
    'wmapN','wmapVolSurgSum','wmapVolSurgProp', ...
    'wmapMean','wmapLCI','wmapUCI','wmapP025','wmapP975', ...
    'petP975_ExtAll','petP975_ExtTemp','wmapP975_ExtAll','wmapP975_ExtTemp','normHippoVol_asym'});

%  Clinical variables and outcome
Side_PTs = Demo_disPTs.Side;  Side_PTs_num = zeros(size(Side_PTs));
Side_PTs_num(strcmp(Side_PTs,'L')) = 1;  Side_PTs_num(strcmp(Side_PTs,'R')) = 2;
Dur_PTs = Demo_disPTs.Duration;
SZFreq_PTs = Demo_disPTs.SZFreq;  SZFreq_PTs_num = zeros(size(SZFreq_PTs));
SZFreq_PTs_num(strcmp(SZFreq_PTs,'Daily')) = 1;   SZFreq_PTs_num(strcmp(SZFreq_PTs,'Weekly')) = 2;
SZFreq_PTs_num(strcmp(SZFreq_PTs,'Monthly')) = 3; SZFreq_PTs_num(strcmp(SZFreq_PTs,'Yearly')) = 4;

Patho_PTs = Demo_disPTs.Pathology;
Patho_bi_num = zeros(size(Patho_PTs));
Patho_bi_num(strcmp(Patho_PTs,'Gliosis')) = 1;
Patho_bi_num(strcmp(Patho_PTs,'HS')) = 2;
Patho_bi_num(strncmp(Patho_PTs,'LITT',4)) = 2;
Patho_bi_num(strncmp(Patho_PTs,'FCD 1',5)) = 3;
Patho_bi_num(strncmp(Patho_PTs,'FCD 3a',6)) = 4;

surgTime = Demo_disPTs.surgTime;

Prog_PTs = Demo_disPTs.x1YProgI;
Prog_bi = repmat({'NSF'},size(Prog_PTs));
Prog_bi(strcmp(Prog_PTs,'ILAE 1')) = {'SF'};
Prog_bi_num = zeros(size(Prog_bi));
Prog_bi_num(strcmp(Prog_bi,'SF')) = 1;  Prog_bi_num(strcmp(Prog_bi,'NSF')) = 0;

dataX = double([Side_PTs_num,Dur_PTs,SZFreq_PTs_num,Patho_bi_num, ...
    surgVolProp,ratio_P975_ExtAll,wmapVolSurgProp,normHippoVol_asym,surgTime]);
dataY = Prog_bi_num;

validMask = surgVolProp >= 0.15;     % NaN (patients not in the new data) -> false
dataX_clean = dataX(validMask,:);
dataY_clean = dataY(validMask);

% Original script removed row 187 of the cleaned data as an outlier. The
% cleaned row order can change with the new data, so check this row.
if size(dataX_clean,1) >= 187
    fprintf('Removing cleaned row 187 as in FIG7B (check this is still the intended patient).\n');
    dataX_clean(187,:) = []; dataY_clean(187) = [];
end
fprintf('FIG7B: %d patients enter the logistic analysis.\n',size(dataX_clean,1));

%  Ridge plot (odds ratios)
x1 = dataX_clean(:,8);
x2 = dataX_clean(:,7);
x3 = dataX_clean(:,5);
x4 = dataX_clean(:,6);
y  = dataY_clean;

data = [x1,x2,x3,x4];
varNames = {'normHippoVol_asym','wmapVolSurgProp','surgVolProp','border'};

group_colors = [cmap_G2(128,:); cmap_G1(128,:); cmap_mean(128,:); blue(128,:)];

binW = 0.01; dotR = 0.020; dotSz = 4; spread = 0.08;

load_paper
figure('Units','centimeters','Position',[2 2 21 29.7]);
ax = subplot(10,7,[1,2,8,9,15,16]);
hold on;
ax.FontSize = 10; ax.LineWidth = 1; ax.XColor = 'k'; ax.YColor = 'k';
ax.YTick = 1:4;
ax.YTickLabel = fliplr(varNames);
ylim([0.3 4.7]);
set(gca,'XScale','log');
xline(1,'--','LineWidth',1,'Color',[0.5 0.5 0.5]);

pred_prob = nan(size(y,1),4);
OR_all = nan(4,1); CI_all = nan(4,2); P_all = nan(4,1);
for v = 1:4
    ypos = 5-v;
    vals_raw = data(:,v);
    vals = zscore(vals_raw);
    col = group_colors(v,:);

    tbl = table(vals,y,'VariableNames',{'x','y'});
    mdl = fitglm(tbl,'y ~ x','Distribution','binomial');
    beta = mdl.Coefficients.Estimate(2);
    se   = mdl.Coefficients.SE(2);
    OR   = exp(beta);
    CI1  = exp(beta-1.96*se);
    CI2  = exp(beta+1.96*se);
    pred_prob(:,v) = mdl.Fitted.Probability;

    zStat = beta/se;
    pVal = 2*(1-normcdf(abs(zStat)));
    OR_all(v) = OR; CI_all(v,:) = [CI1 CI2]; P_all(v) = pVal;

    vals_plot = zscore(vals_raw)*spread+OR;

    [fk,xi] = ksdensity(vals_plot,'NumPoints',256);
    fk = fk(:)/max(fk)*0.35;
    xi = xi(:);
    py = [repmat(ypos,size(xi)); flipud(ypos+fk)];
    px = [xi; flipud(xi)];
    fill(px,py,col,'FaceAlpha',0.65,'EdgeColor',col*0.7);

    edges = (min(vals_plot)-binW):binW:(max(vals_plot)+binW);
    for b = 1:numel(edges)-1
        inBin = vals_plot(vals_plot >= edges(b) & vals_plot < edges(b+1));
        nDots = numel(inBin);
        if nDots == 0, continue; end
        xCenter = mean([edges(b),edges(b+1)]);
        yDots = ypos-dotR*(2*(1:nDots)-1);
        scatter(repmat(xCenter,1,nDots),yDots,dotSz,col,'filled', ...
            'MarkerFaceAlpha',0.80,'MarkerEdgeAlpha',0);
    end

    plot([CI1 CI2],[ypos ypos],'k-','LineWidth',1);
    plot(OR,ypos,'o','MarkerSize',6,'MarkerFaceColor',col,'MarkerEdgeColor','k');

    fprintf('%-25s  OR = %.3f  95%% CI [%.3f, %.3f]  p = %.4f\n', ...
        varNames{v},OR,CI1,CI2,pVal);
end

xlim([0.6 2.3]);
xticks([0.6 0.8 1 1.5 2]);
xticklabels({'0.6','0.8','1','1.5','2'});
xlabel('Odds Ratio');

Results_FIG7B = table(string(varNames(:)),OR_all,CI_all(:,1),CI_all(:,2),P_all, ...
    'VariableNames',{'Variable','OR','CI_low','CI_high','PValue'});
disp(Results_FIG7B);

%  Partial dependence (ICE)
featureNames = {'Side','Dur','SZFreq','Patho','surgVolProp','boundary', ...
    'wmapVolSurgProp','normHippoVol_asym','surgTime'};
X_table = array2table(dataX_clean,'VariableNames',featureNames);
mdl = fitclinear(X_table,dataY_clean,'Learner','logistic');

figure;
plotPartialDependence(mdl,'wmapVolSurgProp',mdl.ClassNames(1),X_table, ...
    'Conditional','absolute');

surgData_clean = [dataX_clean,pred_prob(:,2),dataY_clean];
featureNames2 = [featureNames,'predSF','surgProg'];
disSet_T = array2table(surgData_clean,'VariableNames',featureNames2);

fprintf('\nCompleted eFIG6_12_1 analysis (non-temporal non-insular FCD R-L reference).\n');