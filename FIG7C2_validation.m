clear;clc;

%% 1. Load data, surface, color
load_demo
load_surface
load_color
load('/host/verges/tank/data/mo/demo/02Voxel_featsPET_hmz_flip_asym_wPDB.mat')

Voxel_valPTs_neofeats_hmz_flip_asym_wPDB.MidPET = Voxel_allPTs_neofeats_hmz_flip_asym_wPDB.MidPET([length(Demo_disPTs.ID) + 1: length(Demo_disPTs.ID) + 38], :);
%% 2. Surgical metrices collection
% Inputs 
% PET   : N x V  (continuous)
% Surg  : N x V  (0/1)
% Wabn  : N x V  (0/1)  OR  Wscore : N x V (continuous)
% TL    : 1 x V  (0/1)  (same for all subjects)
N = length(Demo_valPTs.ID);
% preallocate
surgVol = nan(N,1);
surgArea  = nan(N,1);

petMean  = nan(N,1);
petLCI = nan(N,1);
petUCI = nan(N,1);
petP025  = nan(N,1);
petP975  = nan(N,1);
petN     = nan(N,1);

wmapMean  = nan(N,1);
wmapLCI = nan(N,1);
wmapUCI = nan(N,1);
wmapP025  = nan(N,1);
wmapP975  = nan(N,1);

wmapN = nan(N,1);
wmapVolSum    = nan(N,1);
wmapVolSurgSum   = nan(N,1);
wmapVolSurgProp  = nan(N,1);
wmapVolSurgN = nan(N,1);

petP975_ExtAll  = nan(N,1);
petP975_ExtTemp = nan(N,1);
wmapP975_ExtAll = nan(N,1);
wmapP975_ExtTemp = nan(N,1);
ratio_P975_ExtAll = nan(N,1);
ratio_P975_ExtTemp = nan(N,1);

maskAll = double(ismember(economo_parc, allTempLbls));
maskTemp = double(ismember(economo_parc, pureTempLbls));

PATH_valPTs = '/host/verges/tank/data/mo/BIDS_20260316_TLE_BTTH2_BSBH/derivatives/';

for i = 1:length(Demo_valPTs.ID)
    % indi PET
    indi_MidPET_hmz_flip_asym_wPDB = Voxel_valPTs_neofeats_hmz_flip_asym_wPDB.MidPET(i,1:32492)';
    % plot_hemispheres(indi_MidPET_hmz_flip_asym_wPDB, surf_lh, 'views', 'lm'); colormap(flipud(blue))

    % w-score PETmap
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal = indi_MidPET_hmz_flip_asym_wPDB;
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal(~(maskAll(1:end/2, 1))) = 0;
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal(indi_MidPET_hmz_flip_asym_wPDB_ThrVal > -1.96) = 0;
    % plot_hemispheres(indi_MidPET_hmz_flip_asym_wPDB_ThrVal, surf_lh, 'views', 'lm'); colormap(flipud(blue)); 

    % surg 
    File_valPTs_LH_surgMask = strcat(PATH_valPTs,'surgMask/',Demo_valPTs.ID{i},'/ses-post/',Demo_valPTs.ID{i},'_ses-post_hemi-L_space-nativepro_surf-fsLR-32k_label-cavity.func.gii');
    File_valPTs_RH_surgMask = strcat(PATH_valPTs,'surgMask/',Demo_valPTs.ID{i},'/ses-post/',Demo_valPTs.ID{i},'_ses-post_hemi-R_space-nativepro_surf-fsLR-32k_label-cavity.func.gii');
    assume_surg_raw = [gifti(File_valPTs_LH_surgMask).cdata; gifti(File_valPTs_RH_surgMask).cdata];

    if strcmp(Demo_valPTs.Side{i}, 'R')
        assume_surg_temp = [assume_surg_raw(end/2+1:end, 1); assume_surg_raw(1:end/2, 1)];
    else
        assume_surg_temp = assume_surg_raw;
    end
    
    in_mask = ismember(economo_parc, allTempLbls);  
    assume_surg_temp(~in_mask) = 0;               
    assume_surg = assume_surg_temp(1:end/2, 1);

    % subplot(1,4,1); imagesc(assume_surg_raw); subplot(1,4,2); imagesc(assume_surg_temp);
    % subplot(1,4,3); imagesc(in_mask); subplot(1,4,4); imagesc(assume_surg)
    % plot_hemispheres(assume_surg, surf_lh, 'views', 'lm', 'labeltext', Demo_PTs.ID{i}); colormap(white2black);%%%%

    % ---- 1) Resection size within TL ----
    idxTempSurg = (assume_surg) & maskAll(1:end/2,1);
    surgVol(i) = sum(idxTempSurg);
    surgVolProp = surgVol / sum(maskAll(1:end/2,1));
    % ---- 2) PET within resection (TL) + 95% CI ----
    petVals = indi_MidPET_hmz_flip_asym_wPDB(idxTempSurg, 1);
    petVals = petVals(~isnan(petVals));
    petN(i) = numel(petVals);

    if petN(i) >= 2
        mu = mean(petVals);
        sd = std(petVals, 0);
        tcrit = tinv(0.975, petN(i)-1);
        half = tcrit * sd / sqrt(petN(i));

        petMean(i) = mu;
        petLCI(i) = mu - half;    petUCI(i) = mu + half;

        pr = prctile(petVals, [2.5 97.5]);
        petP025(i) = pr(1);       petP975(i) = pr(2);

    elseif petN(i) == 1
        petMean(i) = petVals;
        petLCI(i) = NaN; petUCI(i) = NaN;
        petP025(i) = petVals;      petP975(i) = petVals;
    end

    % ---- 3) Proportion of abnormal w-score resected (TL) ----
    wmapThr = indi_MidPET_hmz_flip_asym_wPDB_ThrVal(:,1) <= -1.96;
    wmapN(i) = sum(wmapThr);
    wmapVals = indi_MidPET_hmz_flip_asym_wPDB_ThrVal(wmapThr,1);
    
    idxWmapTempSurg = wmapThr & idxTempSurg;
    wmapVolSum(i)  = sum(wmapThr);
    wmapVolSurgSum(i) = sum(idxWmapTempSurg);

    if wmapVolSum(i) > 0
        wmapVolSurgProp(i) = wmapVolSurgSum(i) / wmapVolSum(i);
    else
        wmapVolSurgProp(i) = 0; % no abnormal vertices in TL
    end

    % ---- 4) wscore within resection (TL) + 95% CI ----
    wmapVolSurg = indi_MidPET_hmz_flip_asym_wPDB_ThrVal(idxWmapTempSurg,1);
    wmapVolSurg = wmapVolSurg(~isnan(wmapVolSurg));
    wmapVolSurgN(i) = numel(wmapVolSurg);

    if wmapVolSurgN(i) >= 2
        mu = mean(wmapVolSurg);
        sd = std(wmapVolSurg, 0);
        tcrit = tinv(0.975, wmapVolSurgN(i)-1);
        half = tcrit * sd / sqrt(wmapVolSurgN(i));

        wmapMean(i) = mu;
        wmapLCI(i) = mu - half;     wmapUCI(i) = mu + half;

        pr = prctile(wmapVolSurg, [2.5 97.5]);
        wmapP025(i) = pr(1);        wmapP975(i) = pr(2);

    elseif wmapVolSurgN(i) == 1
        wmapMean(i) = wmapVolSurg;
        wmapLCI(i) = NaN; wmapUCI(i) = NaN;
        wmapP025(i) = wmapVolSurg;   wmapP975(i) = wmapVolSurg;
    end

    petExtAllMin = min(indi_MidPET_hmz_flip_asym_wPDB(maskAll == 1,1));
    petExtAllMax = max(indi_MidPET_hmz_flip_asym_wPDB(maskAll == 1,1));

    petExtTempMin = min(indi_MidPET_hmz_flip_asym_wPDB(maskTemp == 1,1));
    petExtTempMax = max(indi_MidPET_hmz_flip_asym_wPDB(maskTemp == 1,1));

    x = petP975(i); 
    [f,xx] = ecdf(indi_MidPET_hmz_flip_asym_wPDB(maskAll == 1, 1));
    p = f(find(xx >= x, 1, 'first'));

    % % figure
    % plot(xx,f, 'b-'); 
    % hold on;
    % plot([x x], [0 p], 'r--')
    % plot([min(xx) x], [p p], 'r--')
    % plot(x,p,'ro')

    petP975_ExtAll(i)  = (petP975(i) - petExtAllMin) / (petExtAllMax - petExtAllMin);
    petP975_ExtTemp(i) = (petP975(i) - petExtTempMin) / (petExtTempMax - petExtTempMin);
    % wmapP975_ExtAll(i) = (wmapP975(i) - petExtAllMin) / (petExtAllMax - petExtAllMin);
    % wmapP975_ExtTemp(i) = (wmapP975(i) - petExtTempMin) / (petExtTempMax - petExtTempMin);
    ratio_P975_ExtAll(i) = petP975(i)  / petExtAllMax;
    ratio_P975_ExtTemp(i) = petP975(i) / petExtTempMax;

end

% Build output table
T = table( (1:N)', surgVol, surgArea, surgVolProp, ...
           petMean, petLCI, petUCI, petP025, petP975, ...
           wmapN, wmapVolSurgN, wmapVolSurgProp, ...
           wmapMean, wmapLCI, wmapUCI, wmapP025, wmapP975, ...
           petP975_ExtAll, petP975_ExtTemp, wmapP975_ExtAll, wmapP975_ExtTemp, ...
           'VariableNames', {'Subj','surgVol','surgArea', 'surgVolProp', ...
                             'petMean','petLCI','petUCI','PET_p2p5','petP975', ...
                             'wmapN','wmapVolSurgSum','wmapVolSurgProp', ...
                             'wmapMean', 'wmapLCI', 'wmapUCI', 'wmapP025', 'wmapP975', ...
                             'petP975_ExtAll', 'petP975_ExtTemp', 'wmapP975_ExtAll', 'wmapP975_ExtTemp'});

%% 4. Sort out data and prepare to analyze
Age_valPTs    = Demo_valPTs.AgeSurg;
Side_valPTs   = Demo_valPTs.Side;     Side_valPTs_num = zeros(size(Side_valPTs));
Side_valPTs_num(strcmp(Side_valPTs, 'L')) = 1;   Side_valPTs_num(strcmp(Side_valPTs, 'R')) = 2;
Dur_valPTs    = Demo_valPTs.Duration;
SZFreq_valPTs = Demo_valPTs.SZFreq;   SZFreq_valPTs_num = zeros(size(SZFreq_valPTs));
SZFreq_valPTs_num(strcmp(SZFreq_valPTs, 'Daily')) = 1;   SZFreq_valPTs_num(strcmp(SZFreq_valPTs, 'Weekly')) = 2;
SZFreq_valPTs_num(strcmp(SZFreq_valPTs, 'Monthly')) = 3; SZFreq_valPTs_num(strcmp(SZFreq_valPTs, 'Yearly')) = 4;

Patho_valPTs = Demo_valPTs.Pathology;
Patho_bi_num = zeros(size(Patho_valPTs));
% Gliosis = 1
Patho_bi_num(strcmp(Patho_valPTs, 'Gliosis')) = 1;
% HS and LITT/NA = 2
Patho_bi_num(strcmp(Patho_valPTs, 'HS')) = 2;
Patho_bi_num(strncmp(Patho_valPTs, 'LITT', 4)) = 2;
% Anything starting with 'FCD' = 3
Patho_bi_num(strncmp(Patho_valPTs, 'FCD 1', 5)) = 3;
Patho_bi_num(strncmp(Patho_valPTs, 'FCD 3a', 6)) = 4;

surgTime   = Demo_valPTs.surgTime;

Prog_valPTs   = Demo_valPTs.x1YProgI; 
Prog_bi    = repmat({'NSF'}, size(Prog_valPTs)); 
Prog_bi(strcmp(Prog_valPTs, 'ILAE 1')) = {'SF'};
Prog_bi_num = zeros(size(Prog_bi));
Prog_bi_num(strcmp(Prog_bi, 'SF')) = 1;   Prog_bi_num(strcmp(Prog_bi, 'NSF')) = 0;

normHippoVol_asym = zeros(length(Demo_valPTs.ID),1);

dataX = double([Side_valPTs_num, Dur_valPTs, SZFreq_valPTs_num, Patho_bi_num, ...
    surgVolProp, ratio_P975_ExtAll, wmapVolSurgProp, normHippoVol_asym, surgTime]);
dataY = Prog_bi_num;

validMask = surgVolProp >= 0.15;
dataX_clean = dataX(validMask, :);
dataY_clean = dataY(validMask);

%% 5. calculate the propSF
y  = dataY_clean;
vals_raw = dataX_clean(:,7);
vals     = zscore(vals_raw);

tbl = table(vals, y, 'VariableNames', {'x','y'});
mdl = fitglm(tbl, 'y ~ x', 'Distribution', 'binomial');
beta  = mdl.Coefficients.Estimate(2);
se    = mdl.Coefficients.SE(2);
OR    = exp(beta);
CI1   = exp(beta - 1.96*se);
CI2   = exp(beta + 1.96*se);
pred_prob = mdl.Fitted.Probability;

%% 6. save as matrix
surgData_clean = [dataX_clean, pred_prob, dataY_clean];
featureNames = {'Side', 'Dur', 'SZFreq', 'Patho', 'surgVolProp','boundary','wmapVolSurgProp', 'normHippoVol_asym', 'surgTime'};
featureNames2 = [featureNames, 'predSF','surgProg'];
valSet_T = array2table(surgData_clean, 'VariableNames', featureNames2);
% save('/host/verges/tank/data/mo/demo/55valSet_T', 'valSet_T')

%% %% for check data, group comparison
% % ── exclude surgVolProp < 0.15 and wmapN == 0 once, applied to all metrics ──
% Prog_valPTs   = Demo_valPTs.x1YProgI; 
% Prog_bi    = repmat({'NSF'}, size(Prog_valPTs)); 
% Prog_bi(strcmp(Prog_valPTs, 'ILAE 1')) = {'SF'};
% Prog_bi_num = zeros(size(Prog_bi));
% Prog_bi_num(strcmp(Prog_bi, 'SF')) = 1;   Prog_bi_num(strcmp(Prog_bi, 'NSF')) = 0;
% 
% 
% validMask = T.surgVolProp >= 0.0 & T.wmapN > 0;
% T_clean   = T(validMask, :);
% Prog_bi_clean = Prog_bi(validMask);
% fprintf('Filter: %d subjects removed, %d retained.\n', ...
%         sum(~validMask), sum(validMask));
% % ────────────────────────────────────────────────────────────────────────────
% 
% sfIdx  = strcmp(Prog_bi_clean, 'SF');
% nsfIdx = strcmp(Prog_bi_clean, 'NSF');
% 
% metrics      = {'surgVolProp', 'petP975_ExtAll', 'petP975_ExtTemp', 'wmapVolSurgProp'};
% metricLabels = {'surgVolProp', 'PET p97.5 (All)', 'PET p97.5 (Temp)', 'Abnormal W-score Resected (%)'};
% figure('Position', [100 100 1200 400]);
% p_temp = zeros(1,4);
% 
% for m = 1:4
%     vals     = T_clean.(metrics{m});
%     sf_vals  = vals(sfIdx);
%     nsf_vals = vals(nsfIdx);
% 
%     % remove NaNs
%     sf_vals  = sf_vals(~isnan(sf_vals));
%     nsf_vals = nsf_vals(~isnan(nsf_vals));
% 
%     % t-test
%     [h, p_val, ~, stats] = ttest2(sf_vals, nsf_vals);
%     subplot(1, 4, m);
% 
%     % boxplot
%     grp_data = [sf_vals;  nsf_vals];
%     grp_lbls = [repmat({'SF'},  numel(sf_vals),  1); ...
%                 repmat({'NSF'}, numel(nsf_vals), 1)];
%     boxplot(grp_data, grp_lbls, 'Colors', 'kb', 'Symbol', 'o');
% 
%     % overlay individual points with jitter
%     hold on;
%     jitter = 0.07;
%     scatter(ones(size(sf_vals))   + (rand(size(sf_vals)) -0.5)*jitter,  sf_vals,  30, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
%     scatter(2*ones(size(nsf_vals))+ (rand(size(nsf_vals))-0.5)*jitter, nsf_vals, 30, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
% 
%     % p-value annotation
%     ymax   = max([sf_vals; nsf_vals]);
%     ymin   = min([sf_vals; nsf_vals]);
%     yrange = ymax - ymin;
%     line([1 2], [ymax + 0.08*yrange, ymax + 0.08*yrange], 'Color', 'k');
%     if p_val < 0.001
%         sig_str = 'p < 0.001';
%     elseif p_val < 0.01
%         sig_str = sprintf('p = %.3f **', p_val);
%     elseif p_val < 0.05
%         sig_str = sprintf('p = %.3f *',  p_val);
%     else
%         sig_str = sprintf('p = %.3f ns', p_val);
%     end
%     text(1.5, ymax + 0.12*yrange, sig_str, 'HorizontalAlignment', 'center', 'FontSize', 10);
%     ylabel(metricLabels{m});
% 
%     fprintf('--- %s ---\n', metrics{m});
%     fprintf('  SF:  mean=%.3f, SD=%.3f, n=%d\n', mean(sf_vals),  std(sf_vals),  numel(sf_vals));
%     fprintf('  NSF: mean=%.3f, SD=%.3f, n=%d\n', mean(nsf_vals), std(nsf_vals), numel(nsf_vals));
%     fprintf('  t(%.0f)=%.3f, p=%.4f\n\n', stats.df, stats.tstat, p_val);
%     p_temp(1,m) = p_val;
% end
% 
% [h_prog, crit_p_prog, adj_p_prog] = fdr_bh(p_temp, 0.05, 'pdep', 'yes');
% display(h_prog);

%% 7. disSet and valSet
%  7.1 Load and sort both datasets
load('/host/verges/tank/data/mo/demo/55disSet_T.mat')
dis_T = sortrows(disSet_T, 'predSF', 'ascend');
load('/host/verges/tank/data/mo/demo/55valSet_T.mat')
val_T = sortrows(valSet_T, 'predSF', 'ascend');
N_dis = height(dis_T);
N_val = height(val_T);
new_T = [dis_T; val_T];
N     = height(new_T);

% Virtual X layout: Total virtual range: [0, 100]
x_total = 100;
dis_start = 0;   dis_end = 60;
gap_start = 60;  gap_end = 65;
val_start = 65;  val_end = 100;
div_x     = gap_start;

x_dis  = linspace(dis_start, dis_end, N_dis)';
x_val  = linspace(val_start, val_end, N_val)';
x_virt = [x_dis; x_val];

%  7.2 Heatmap data — pick ONE or MULTIPLE variables to display
% >>>  CHANGE THIS LINE to switch variables  <
showVars = {'wmapVolSurgProp'};   % add/remove names freely
% options: 'normHippoVol_asym' | 'wmapVolSurgProp' | 'surgVolProp' | 'boundary'

geneNames = {'normHippoVol_asym','wmapVolSurgProp','surgVolProp','boundary'};
G  = [new_T.normHippoVol_asym, new_T.wmapVolSurgProp, ...
      new_T.surgVolProp,        new_T.boundary];

% Automatically find columns for selected variables and z-score each
varIdx = find(ismember(geneNames, showVars));
varIdx = varIdx(ismember(geneNames(varIdx), showVars));  % preserve showVars order
% re-order to match showVars order (not geneNames order)
[~, orderIdx] = ismember(showVars, geneNames);
G_sel = G(:, orderIdx);                        % [N x nVars]

nVars = numel(showVars);
Gz = zeros(N, nVars);
for j = 1:nVars
    mu = mean(G_sel(:,j), 'omitnan');
    sd = std( G_sel(:,j), 0, 'omitnan');
    if sd == 0, sd = 1; end
    Gz(:,j) = (G_sel(:,j) - mu) ./ sd;        % [N x nVars]
end

aliveIdx = new_T.surgProg == 0;
deadIdx  = new_T.surgProg == 1;

fig = figure('Color','w', 'Position',[100 80 760 820]);
tl  = tiledlayout(fig, 4, 1, 'TileSpacing','compact', 'Padding','compact');
xlim_range = [dis_start, val_end];

%  7.3 Risk score
ax1 = nexttile(tl,1);
hold(ax1,'on');

predSF_vals = new_T.predSF;
predSF_vals(N_dis+1:end) = smoothdata(predSF_vals(N_dis+1:end), 'movmean', 3);
pmin = min(predSF_vals,[],'omitnan');
pmax = max(predSF_vals,[],'omitnan');

scatter(ax1, x_virt, predSF_vals, 16, predSF_vals, 'filled');
colormap(ax1, cmap_gray2red);
clim(ax1, [pmin pmax]);

text(ax1, mean([dis_start dis_end]), max(predSF_vals)*1.02, 'DisSet', ...
    'HorizontalAlignment','center', 'FontSize',11, 'FontWeight','bold');
text(ax1, mean([val_start val_end]), max(predSF_vals)*1.02, 'ValSet',  ...
    'HorizontalAlignment','center', 'FontSize',11, 'FontWeight','bold');
ylabel(ax1, 'Risk Score', 'FontSize', 16);
set(ax1, 'XTick',[], 'Box','off', 'LineWidth',1, 'FontSize',12);
ax1.XColor = 'none';  ax1.TickDir = 'out';
xlim(ax1, xlim_range);
yl = ylim(ax1);
ylim(ax1, [yl(1)-0.1*range(yl), yl(2)+0.10*range(yl)]);

%  7.4 Heatmap — one row per selected variable
ax2 = nexttile(tl,2);
hold(ax2,'on');

imagesc(ax2, [dis_start dis_end], [1 nVars], Gz(1:N_dis,:)');
imagesc(ax2, [val_start val_end], [1 nVars], Gz(N_dis+1:end,:)');
colormap(ax2, cmap_gray2red);
clim(ax2, [0 3]);

xline(ax2, div_x, ':k', 'LineWidth', 1);
set(ax2, 'YTick', 1:nVars, 'YTickLabel', showVars, ...  % labels update automatically
    'XTick',[], 'Box','off', 'LineWidth',1, 'FontSize',13, 'TickDir','out');
xlim(ax2, xlim_range);
ylim(ax2, [0.5 nVars+0.5]);                             % adjusts automatically to nVars

cb = colorbar(ax2, 'eastoutside');
cb.Ticks = 0:1:3;
cb.Label.String = 'Expression';
cb.FontSize = 11;

%  7.5 Survival time
ax3 = nexttile(tl,3);
hold(ax3,'on');

scatter(ax3, x_virt(aliveIdx), new_T.surgTime(aliveIdx), 16, gray(16,:), 'filled', ...
    'MarkerEdgeColor','k', 'LineWidth',0.5);
scatter(ax3, x_virt(deadIdx),  new_T.surgTime(deadIdx),  16, red_mean,   'filled', ...
    'MarkerEdgeColor','k', 'LineWidth',0.5);
xline(ax3, div_x, ':k', 'LineWidth', 1.5);
ylabel(ax3, 'Time since surgery', 'FontSize', 16);
set(ax3, 'XTick',[], 'Box','off', 'LineWidth',1.2, 'FontSize',12);
ax3.XColor = 'none';  ax3.TickDir = 'out';
xlim(ax3, xlim_range);

%% 8. example
clear;clc;
load_surface
load_demo
load_color

load('/host/verges/tank/data/mo/demo/54Voxel_valfeats_hmz_flip_asym_wPDB.mat')
PATH_valPTs = '/host/verges/tank/data/mo/BIDS_20260316_TLE_MICs_BTTH2/derivatives/';

maskTemp = double(ismember(economo_parc, pureTempLbls));

%% good example is PX2501016 (i = 22); bad is PX2501053 (i = 4)
for i = 22
    % indi PET
    indi_MidPET_hmz_flip_asym_wPDB = Voxel_valPTs_neofeats_hmz_flip_asym_wPDB.MidPET(i,1:32492)';
    % plot_hemispheres(indi_MidPET_hmz_flip_asym_wPDB, surf_lh, 'views', 'lm'); colormap(flipud(blue))

    % w-score PETmap
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal = indi_MidPET_hmz_flip_asym_wPDB;
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal(~(maskTemp(1:end/2, 1))) = 0;
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal(indi_MidPET_hmz_flip_asym_wPDB_ThrVal > -1.96) = 0;
    % plot_hemispheres(indi_MidPET_hmz_flip_asym_wPDB_ThrVal, surf_lh, 'views', 'lm'); colormap(flipud(blue)); 

    % surg 
    File_PTs_LH_surgMask = strcat(PATH_valPTs,'surgMask/',Demo_valPTs.ID{i},'/ses-post/',Demo_valPTs.ID{i},'_ses-post_hemi-L_space-nativepro_surf-fsLR-32k_label-cavity.func.gii');
    File_PTs_RH_surgMask = strcat(PATH_valPTs,'surgMask/',Demo_valPTs.ID{i},'/ses-post/',Demo_valPTs.ID{i},'_ses-post_hemi-R_space-nativepro_surf-fsLR-32k_label-cavity.func.gii');
    assume_surg_raw = [gifti(File_PTs_LH_surgMask).cdata; gifti(File_PTs_RH_surgMask).cdata];

    if strcmp(Demo_valPTs.Side{i}, 'R')
        assume_surg_temp = [assume_surg_raw(end/2+1:end, 1); assume_surg_raw(1:end/2, 1)];
    else
        assume_surg_temp = assume_surg_raw;
    end
    
    in_mask = ismember(economo_parc, allTempLbls);  
    assume_surg_temp(~in_mask) = 0;               
    assume_surg = assume_surg_temp(1:end/2, 1);
    
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal(indi_MidPET_hmz_flip_asym_wPDB_ThrVal ~= 0) = 1;
    plot_hemispheres([indi_MidPET_hmz_flip_asym_wPDB, indi_MidPET_hmz_flip_asym_wPDB_ThrVal, assume_surg], surf_lh, ...
        'labeltext', {Demo_valPTs.ID(i), Demo_valPTs.x1YProgI(i), string(valSet_T.wmapVolSurgProp(i))}, 'views', 'lm'); colormap(flipud(blue)); 
    % plot_hemispheres(indi_MidPET_hmz_flip_asym_wPDB_ThrVal, surf_lh, 'views', 'lm'); colormap([1 1 1; red_mean]);
    % plot_hemispheres(assume_surg, surf_lh, 'views', 'lm'); colormap(cmap_purple); 
    % plot_hemispheres(assume_surg, surf_lh, 'views', 'lm'); colormap(white2black);
    % plot_hemispheres(indi_MidPET_hmz_flip_asym_wPDB, surf_lh, 'views', 'lm'); colormap(flipud(blue)); 
end