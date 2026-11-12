clear;clc;

%% 1. Load data, surface, color
load_demo
load_surface
load_color

%% 2. Surgical metrices collection
% Inputs 
% PET   : N x V  (continuous)
% Surg  : N x V  (0/1)
% Wabn  : N x V  (0/1)  OR  Wscore : N x V (continuous)
% TL    : 1 x V  (0/1)  (same for all subjects)
N = length(Demo_disPTs.ID);
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

%  2.2 Individual participant Voxel_MidPET_hmz_flip_asym_wPDB
load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat')
PATH_PTs = '/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/derivatives/';

for i = 1:length(Demo_disPTs.ID)
    % indi PET
    indi_MidPET_hmz_flip_asym_wPDB = Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET(i,1:32492)';
    % plot_hemispheres(indi_MidPET_hmz_flip_asym_wPDB, surf_lh, 'views', 'lm'); colormap(flipud(blue))

    % w-score PETmap
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal = indi_MidPET_hmz_flip_asym_wPDB;
    indi_MidPET_hmz_flip_asym_wPDB_ThrVal(~(maskAll(1:end/2, 1))) = 0;
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

    % subplot(1,4,1); imagesc(assume_surg_raw); subplot(1,4,2); imagesc(assume_surg_temp);
    % subplot(1,4,3); imagesc(in_mask); subplot(1,4,4); imagesc(assume_surg)
    % plot_hemispheres(assume_surg, surf_lh, 'views', 'lm', 'labeltext', Demo_disPTs.ID{i}); colormap(white2black);%%%%

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

%% 3. Read hippocampal volume data and extract normalized volumes
hipp_data = readtable('/host/verges/tank/data/mo/demo/52all_hippocampal_volumes_normalized.csv');

n_patients = height(Demo_disPTs);

% Preallocate
normHippoVol_IPS  = nan(n_patients, 1);
normHippoVol_CON  = nan(n_patients, 1);
normHippoVol_asym = nan(n_patients, 1);

for i = 1:n_patients
    subj = Demo_disPTs.ID{i};
    side = Demo_disPTs.Side{i};  % 'L' or 'R'

    % Contralateral side
    if strcmp(side, 'L')
        contra_side = 'R';
    else
        contra_side = 'L';
    end

    subj_rows = strcmp(hipp_data.subject, subj);

    if sum(subj_rows) == 0
        warning('Subject %s not found in hippocampal volume file', subj);
        continue
    end

    % Ipsilateral
    ipsi_row  = subj_rows & strcmp(hipp_data.hemi, side);
    % Contralateral
    contra_row = subj_rows & strcmp(hipp_data.hemi, contra_side);

    if sum(ipsi_row) == 0 || sum(contra_row) == 0
        warning('Missing hemisphere for subject %s', subj);
        continue
    end

    normHippoVol_IPS(i)  = hipp_data.normalized_hipp_vol(ipsi_row);
    normHippoVol_CON(i)  = hipp_data.normalized_hipp_vol(contra_row);

    % Asymmetry index (same formula as your parcellation code)
    normHippoVol_asym(i) = (normHippoVol_IPS(i) - normHippoVol_CON(i)) / ...
                           (0.5 * (abs(normHippoVol_IPS(i)) + abs(normHippoVol_CON(i))));
end

% Build output table
T = table( (1:N)', surgVol, surgArea, surgVolProp, ...
           petMean, petLCI, petUCI, petP025, petP975, ...
           wmapN, wmapVolSurgN, wmapVolSurgProp, ...
           wmapMean, wmapLCI, wmapUCI, wmapP025, wmapP975, ...
           petP975_ExtAll, petP975_ExtTemp, wmapP975_ExtAll, wmapP975_ExtTemp, normHippoVol_asym, ...
           'VariableNames', {'Subj','surgVol','surgArea', 'surgVolProp', ...
                             'petMean','petLCI','petUCI','PET_p2p5','petP975', ...
                             'wmapN','wmapVolSurgSum','wmapVolSurgProp', ...
                             'wmapMean', 'wmapLCI', 'wmapUCI', 'wmapP025', 'wmapP975', ...
                             'petP975_ExtAll', 'petP975_ExtTemp', 'wmapP975_ExtAll', 'wmapP975_ExtTemp', 'normHippoVol_asym'});

%% 4. Sort out data and prepare to analyze
Age_PTs    = Demo_disPTs.AgeSurg;
Side_PTs   = Demo_disPTs.Side;     Side_PTs_num = zeros(size(Side_PTs));
Side_PTs_num(strcmp(Side_PTs, 'L')) = 1;   Side_PTs_num(strcmp(Side_PTs, 'R')) = 2;
Dur_PTs    = Demo_disPTs.Duration;
SZFreq_PTs = Demo_disPTs.SZFreq;   SZFreq_PTs_num = zeros(size(SZFreq_PTs));
SZFreq_PTs_num(strcmp(SZFreq_PTs, 'Daily')) = 1;   SZFreq_PTs_num(strcmp(SZFreq_PTs, 'Weekly')) = 2;
SZFreq_PTs_num(strcmp(SZFreq_PTs, 'Monthly')) = 3; SZFreq_PTs_num(strcmp(SZFreq_PTs, 'Yearly')) = 4;

% Patho_PTs = Demo_disPTs.Pathology;
% Patho_bi  = repmat({'EpiPatho'}, size(Patho_PTs)); 
% Patho_bi(strcmp(Patho_PTs, 'Gliosis')) = {'Gliosis'};
% Patho_bi_num = zeros(size(Patho_bi));
% Patho_bi_num(strcmp(Patho_bi, 'Gliosis')) = 1;   Patho_bi_num(strcmp(Patho_bi, 'EpiPatho')) = 2;
Patho_PTs = Demo_disPTs.Pathology;
Patho_bi_num = zeros(size(Patho_PTs));
% Gliosis = 1
Patho_bi_num(strcmp(Patho_PTs, 'Gliosis')) = 1;
% HS and LITT/NA = 2
Patho_bi_num(strcmp(Patho_PTs, 'HS')) = 2;
Patho_bi_num(strncmp(Patho_PTs, 'LITT', 4)) = 2;
% Anything starting with 'FCD' = 3
Patho_bi_num(strncmp(Patho_PTs, 'FCD 1', 5)) = 3;
Patho_bi_num(strncmp(Patho_PTs, 'FCD 3a', 6)) = 4;

surgTime   = Demo_disPTs.surgTime;

Prog_PTs   = Demo_disPTs.x1YProgI; 
Prog_bi    = repmat({'NSF'}, size(Prog_PTs)); 
Prog_bi(strcmp(Prog_PTs, 'ILAE 1')) = {'SF'};
Prog_bi_num = zeros(size(Prog_bi));
Prog_bi_num(strcmp(Prog_bi, 'SF')) = 1;   Prog_bi_num(strcmp(Prog_bi, 'NSF')) = 0;

dataX = double([Side_PTs_num, Dur_PTs, SZFreq_PTs_num, Patho_bi_num, surgVolProp, ratio_P975_ExtAll, wmapVolSurgProp, normHippoVol_asym, surgTime]);
dataY = Prog_bi_num;

validMask = surgVolProp >= 0.15;
dataX_clean = dataX(validMask, :);
dataY_clean = dataY(validMask);

dataX_clean(187,:) = []; dataY_clean(187)  = [];

%% 5. Ridge plot
x1 = dataX_clean(:,8);
x2 = dataX_clean(:,7);
x3 = dataX_clean(:,5);
x4 = dataX_clean(:,6);
y  = dataY_clean;

data  = [x1, x2, x3, x4];
varNames = {'normHippoVol_asym','wmapVolSurgProp','surgVolProp','border'};

group_colors = [cmap_G2(128,:);
                cmap_G1(128,:);
                cmap_mean(128,:);
                blue(128,:)];

binW  = 0.01;
dotR  = 0.020;
dotSz = 4;
spread = 0.08; % half-width of cloud in log-space units

load_paper
ax = subplot(10, 7, [1,2,8,9,15,16]);
hold on;
ax.FontSize  = 10;
ax.LineWidth = 1;
ax.XColor    = 'k';
ax.YColor    = 'k';
ax.YTick    = 1:4;
ax.YTickLabel = fliplr(varNames);
ylim([0.3 4.7]);
set(gca,'XScale','log');
xline(1,'--','LineWidth',1,'Color',[0.5 0.5 0.5]);

pred_prob = nan(size(y,1),4);
for v = 1:4
    ypos = 5 - v;
    vals_raw = data(:, v);
    vals     = zscore(vals_raw);
    col      = group_colors(v, :);

    % --- logistic regression first to get OR ---
    tbl = table(vals, y, 'VariableNames', {'x','y'});
    mdl = fitglm(tbl, 'y ~ x', 'Distribution', 'binomial');
    beta  = mdl.Coefficients.Estimate(2);
    se    = mdl.Coefficients.SE(2);
    OR    = exp(beta);
    CI1   = exp(beta - 1.96*se);
    CI2   = exp(beta + 1.96*se);
    pred_prob(:,v) = mdl.Fitted.Probability;
    
    % Compute p-value manually from z-score
    zStat = beta / se;
    pVal  = 2 * (1 - normcdf(abs(zStat)));

    % --- center cloud on OR (in log scale) ---
    % zscore vals, then scale spread, then shift center to OR
    vals_plot = zscore(vals_raw) * spread + OR;

    % --- Half violin ---
    [fk, xi] = ksdensity(vals_plot, 'NumPoints', 256);
    fk = fk(:) / max(fk) * 0.35;
    xi = xi(:);
    py = [repmat(ypos, size(xi)); flipud(ypos + fk)];
    px = [xi; flipud(xi)];
    fill(px, py, col, 'FaceAlpha', 0.65, 'EdgeColor', col*0.7);

    % --- Stacked dots ---
    edges = (min(vals_plot)-binW) : binW : (max(vals_plot)+binW);
    for b = 1:numel(edges)-1
        inBin = vals_plot(vals_plot >= edges(b) & vals_plot < edges(b+1));
        nDots = numel(inBin);
        if nDots == 0, continue; end
        xCenter = mean([edges(b), edges(b+1)]);
        yDots   = ypos - dotR*(2*(1:nDots) - 1);
        scatter(repmat(xCenter,1,nDots), yDots, dotSz, col, 'filled', 'MarkerFaceAlpha', 0.80, 'MarkerEdgeAlpha', 0);
    end

    % --- OR + CI ---
    by = ypos + 0.0;
    plot([CI1 CI2], [by by], 'k-', 'LineWidth', 1);
    plot(OR, by, 'o', 'MarkerSize', 6, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k');

    fprintf('%-25s  OR = %.3f  95%% CI [%.3f, %.3f]  p = %.4f\n', varNames{v}, OR, CI1, CI2, pVal);
end

xlim([0.6 2.3]);
xticks([0.6 0.8 1 1.5 2]);
xticklabels({'0.6','0.8','1','1.5','2'});
xlabel('Odds Ratio');

%% 6. PartialDependence plot
featureNames = {'Side', 'Dur', 'SZFreq', 'Patho', 'surgVolProp','boundary','wmapVolSurgProp', 'normHippoVol_asym', 'surgTime'};
X_table = array2table(dataX_clean, 'VariableNames', featureNames);
mdl = fitclinear(X_table, dataY_clean, 'Learner', 'logistic');

% [predictedLabels, scores] = predict(mdl, X_table);
% probY1 = scores(:, 2);  % Probability of class

% % --- PDP ---
% figure;
% plotPartialDependence(mdl, 'wmapVolSurgProp', mdl.ClassNames(1), X_table);
% title('PDP: wmapVolSurgProp → surgProp');

% --- ICE ---
% load_paper
% ax = subplot(7,10,[1 2 8 9]);
figure;
plotPartialDependence(mdl, 'wmapVolSurgProp', mdl.ClassNames(1), X_table, ...
    'Conditional', 'absolute');

%% 7. save
% dataX_clean(:,10) = dataX_clean(:,5)*10;
% dataX_clean(:,11) = dataX_clean(:,7)*10;
surgData_clean = [dataX_clean, pred_prob(:,2), dataY_clean];

featureNames2 = [featureNames, 'predSF','surgProg'];
disSet_T = array2table(surgData_clean, 'VariableNames', featureNames2);
% writetable(disSet_T, '/host/verges/tank/data/mo/demo/52surgData_clean.csv');
% save('/host/verges/tank/data/mo/demo/55disSet_T', 'disSet_T')
