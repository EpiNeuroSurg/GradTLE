clear; clc;

%% 1. Load data
load_demo
load_surface
load_color

load('/host/verges/tank/data/mo/demo/02Voxel_featsPET_hmz_flip_asym_wPDB.mat')

%% 2. Combine variables: ATL = discovery + validation; SAH/LITT = sensitivity

Age_allPTs = [
    Demo_disPTs.AgeSurg;
    Demo_valPTs.AgeSurg;
    Demo_senPTs.AgeSurg
];

Sex_allPTs = [
    Demo_disPTs.Sex;
    Demo_valPTs.Sex;
    Demo_senPTs.Sex
];

Dur_allPTs = [
    Demo_disPTs.Duration;
    Demo_valPTs.Duration;
    Demo_senPTs.Duration
];

SZFreq_allPTs = [
    Demo_disPTs.SZFreq;
    Demo_valPTs.SZFreq;
    Demo_senPTs.SZFreq
];

SurgGrp_allPTs = [
    repmat({'ATL'}, height(Demo_disPTs), 1);
    repmat({'ATL'}, height(Demo_valPTs), 1);
    repmat({'SAH_LITT'}, height(Demo_senPTs), 1)
];

%% 3. Clean seizure frequency

SZFreq_allPTs_num = zeros(size(SZFreq_allPTs));

SZFreq_allPTs_num(strcmpi(SZFreq_allPTs,'Daily'))   = 1;
SZFreq_allPTs_num(strcmpi(SZFreq_allPTs,'Weekly'))  = 2;
SZFreq_allPTs_num(strcmpi(SZFreq_allPTs,'Monthly')) = 3;
SZFreq_allPTs_num(strcmpi(SZFreq_allPTs,'Yearly'))  = 4;

%% 4. Check numbers

fprintf('\n===== Surgical approach cohorts =====\n');
fprintf('All patients: n = %d\n', numel(SurgGrp_allPTs));
fprintf('ATL: n = %d\n', sum(strcmp(SurgGrp_allPTs,'ATL')));
fprintf('SAH/LITT: n = %d\n', sum(strcmp(SurgGrp_allPTs,'SAH_LITT')));

%% 5. Load PET feature

MidPET_allPTs = Voxel_allPTs_neofeats_hmz_flip_asym_wPDB.MidPET;

fprintf('\nPET matrix: %d subjects x %d vertices\n', ...
    size(MidPET_allPTs,1), size(MidPET_allPTs,2));

if size(MidPET_allPTs,1) ~= numel(SurgGrp_allPTs)
    error('Subject number mismatch: PET matrix has %d rows, demographic variables have %d subjects.', ...
        size(MidPET_allPTs,1), numel(SurgGrp_allPTs));
end

%% 6. Define ipsilateral temporal cortex mask using allTempLbls from load_surface

nV = size(MidPET_allPTs,2);
nHalf = nV / 2;

tempMask = false(1,nV);
tempMask(1:nHalf) = ismember(economo_parc(1:nHalf), allTempLbls);

fprintf('\nTemporal labels used: %s\n', mat2str(allTempLbls));
fprintf('Temporal cortex vertices used for FDR: %d\n', sum(tempMask));

%% 7. Surface-wise GLM: ATL vs SAH/LITT

Grp_Cov    = term(SurgGrp_allPTs);
Age_Cov    = term(Age_allPTs);
Sex_Cov    = term(Sex_allPTs);
Dur_Cov    = term(Dur_allPTs);
SZFreq_Cov = term(SZFreq_allPTs_num);

Model = 1 + Grp_Cov + Age_Cov + Sex_Cov + Dur_Cov + SZFreq_Cov;

slm_surg = SurfStatLinMod(MidPET_allPTs, Model);
slm_surg = SurfStatT(slm_surg, Grp_Cov.ATL - Grp_Cov.SAH_LITT);

p_surg = 2 * tcdf(-abs(slm_surg.t), slm_surg.df);
D_surg = slm_surg.ef;

%% 8. FDR correction within temporal cortex only

p_temp = p_surg(tempMask);

[h_temp, crit_p_temp, adj_p_temp] = fdr_bh(p_temp, 0.015, 'pdep', 'yes');

h_surg_temp = false(size(p_surg));
h_surg_temp(tempMask) = h_temp;

idxATL = strcmp(SurgGrp_allPTs,'ATL');
idxSAHLITT = strcmp(SurgGrp_allPTs,'SAH_LITT');

mean_ATL_temp = mean(MidPET_allPTs(idxATL,tempMask), 'all', 'omitnan');
mean_SAH_temp = mean(MidPET_allPTs(idxSAHLITT,tempMask), 'all', 'omitnan');

fprintf('\n===== ATL vs SAH/LITT, temporal cortex only =====\n');
fprintf('Mean HypoPET ATL = %.4f\n', mean_ATL_temp);
fprintf('Mean HypoPET SAH/LITT = %.4f\n', mean_SAH_temp);
fprintf('ATL - SAH/LITT = %.4f\n', mean_ATL_temp - mean_SAH_temp);
fprintf('Minimum temporal uncorrected P = %.6f\n', min(p_temp));
fprintf('Temporal FDR critical P = %.6f\n', crit_p_temp);
fprintf('Number of temporal FDR-significant vertices = %d\n', sum(h_surg_temp));

%% 9. Plot maps

% 9.1 ATL mean HypoPET
PET_ATL = mean(MidPET_allPTs(idxATL,:),1);
PET_ATL(~tempMask) = 0;

plot_hemispheres(PET_ATL(1,1:nHalf)',{surf_lh});
colormap(flipud(blue));
set(findall(gcf,'type','axes'),'CLim',[-6 0]);
title('ATL');

% 9.2 SAH/LITT mean HypoPET
PET_SAHLITT = mean(MidPET_allPTs(idxSAHLITT,:),1);
PET_SAHLITT(~tempMask) = 0;

plot_hemispheres(PET_SAHLITT(1,1:nHalf)',{surf_lh});
colormap(flipud(blue));
set(findall(gcf,'type','axes'),'CLim',[-6 0]);
title('SAH/LITT');

% 9.3 Raw Cohen''s d, temporal cortex only
% D_surg_temp = D_surg;
% D_surg_temp(~tempMask) = 0;
% 
% plot_hemispheres(D_surg_temp(1,1:nHalf)',{surf_lh});
% colormap(flipud(blue));
% set(findall(gcf,'type','axes'),'CLim',[min(D_surg_temp(1,1:nHalf)') 0]);
% title('ATL vs SAH/LITT, Cohen''s d');

D_surg_temp = D_surg;
D_surg_temp(~tempMask) = 0;

plot_hemispheres(D_surg_temp(1,1:nHalf)',{surf_lh});
colormap(flipud(blue));
set(findall(gcf,'type','axes'),'CLim',[-6 0]);
title('ATL vs SAH/LITT, Cohen''s d');

% 9.4 FDR-corrected Cohen''s d
SigMap = zeros(size(D_surg));
SigMap(h_surg_temp) = -1;

plot_hemispheres(SigMap(1,1:nHalf)',{surf_lh});

colormap(flipud(gray));
set(findall(gcf,'type','axes'),'CLim',[-1 0]);

title('ATL vs SAH/LITT, temporal PFDR < 0.05');