clear; clc;

%% 1. Load data
load_demo
load_surface
load_color

load('/host/verges/tank/data/mo/demo/02Voxel_featsPET_hmz_flip_asym_wPDB.mat')

%% 2. Sensitivity cohort only: last 56 patients

nDis = height(Demo_disPTs);
nVal = height(Demo_valPTs);
nSen = height(Demo_senPTs);

sen_idx = (nDis + nVal + 1):(nDis + nVal + nSen);

MidPET_senPTs = Voxel_allPTs_neofeats_hmz_flip_asym_wPDB.MidPET(sen_idx,:);

fprintf('\n===== Sensitivity cohort only =====\n');
fprintf('Sensitivity cohort: n = %d\n', nSen);
fprintf('PET matrix: %d subjects x %d vertices\n', size(MidPET_senPTs,1), size(MidPET_senPTs,2));

%% 3. Demographic variables

Age_senPTs = Demo_senPTs.AgeSurg;
Sex_senPTs = Demo_senPTs.Sex;
Dur_senPTs = Demo_senPTs.Duration;
SZFreq_senPTs = Demo_senPTs.SZFreq;

SZFreq_senPTs_num = zeros(size(SZFreq_senPTs));
SZFreq_senPTs_num(strcmpi(SZFreq_senPTs,'Daily'))   = 1;
SZFreq_senPTs_num(strcmpi(SZFreq_senPTs,'Weekly'))  = 2;
SZFreq_senPTs_num(strcmpi(SZFreq_senPTs,'Monthly')) = 3;
SZFreq_senPTs_num(strcmpi(SZFreq_senPTs,'Yearly'))  = 4;

%% 4. Define SF vs NSF based on Step5A
% 
% Prog_senPTs = Demo_senPTs.x1YProgE;
% 
% Prog_bi = repmat({'NSF'}, size(Prog_senPTs));
% Prog_bi(strcmp(Prog_senPTs, 'Engel 1')) = {'SF'};
% 
% fprintf('SF: n = %d\n', sum(strcmp(Prog_bi,'SF')));
% fprintf('NSF: n = %d\n', sum(strcmp(Prog_bi,'NSF')));
Prog_senPTs = Demo_senPTs.x1YProgE;

Prog_bi = repmat({'NSF'}, size(Prog_senPTs));

Prog_bi(strncmp(Prog_senPTs,'Engel 1',7)) = {'SF'};

fprintf('SF: n = %d\n', sum(strcmp(Prog_bi,'SF')));
fprintf('NSF: n = %d\n', sum(strcmp(Prog_bi,'NSF')));

%% 5. Temporal cortex mask using allTempLbls

nV = size(MidPET_senPTs,2);
nHalf = nV/2;

tempMask = false(1,nV);
tempMask(1:nHalf) = ismember(economo_parc(1:nHalf), allTempLbls);

fprintf('\nTemporal labels used: %s\n', mat2str(allTempLbls));
fprintf('Temporal cortex vertices used for FDR: %d\n', sum(tempMask));

%% 6. Surface-wise GLM: SF vs NSF in sensitivity cohort

Grp_Cov    = term(Prog_bi);
Age_Cov    = term(Age_senPTs);
Sex_Cov    = term(Sex_senPTs);
Dur_Cov    = term(Dur_senPTs);
SZFreq_Cov = term(SZFreq_senPTs_num);

Model = 1 + Grp_Cov + Age_Cov + Sex_Cov + Dur_Cov + SZFreq_Cov;

slm_sen = SurfStatLinMod(MidPET_senPTs, Model);
slm_sen = SurfStatT(slm_sen, Grp_Cov.SF - Grp_Cov.NSF);

p_sen = tcdf(slm_sen.t, slm_sen.df);
D_sen = slm_sen.ef;

%% 7. FDR correction within allTempLbls only

p_temp = p_sen(tempMask);

[h_temp, crit_p_temp, adj_p_temp] = fdr_bh(p_temp, 0.05, 'pdep', 'yes');

h_sen_temp = false(size(p_sen));
h_sen_temp(tempMask) = h_temp;

fprintf('\n===== SF vs NSF, sensitivity cohort, allTempLbls only =====\n');
fprintf('Minimum temporal uncorrected P = %.6f\n', min(p_temp));
fprintf('Temporal FDR critical P = %.6f\n', crit_p_temp);
fprintf('Number of temporal FDR-significant vertices = %d\n', sum(h_sen_temp));

%% 8. Plot maps

idxSF = strcmp(Prog_bi,'SF');
idxNSF = strcmp(Prog_bi,'NSF');

% 8.1 SF mean HypoPET
PET_SF = mean(MidPET_senPTs(idxSF,:),1);
PET_SF(~tempMask) = 0;

plot_hemispheres(PET_SF(1,1:nHalf)',{surf_lh});
colormap(flipud(blue));
set(findall(gcf,'type','axes'),'CLim',[-6 0]);
title('Sensitivity cohort: SF');

% 8.2 NSF mean HypoPET
PET_NSF = mean(MidPET_senPTs(idxNSF,:),1);
PET_NSF(~tempMask) = 0;

plot_hemispheres(PET_NSF(1,1:nHalf)',{surf_lh});
colormap(flipud(blue));
set(findall(gcf,'type','axes'),'CLim',[-6 0]);
title('Sensitivity cohort: NSF');

% 8.3 Raw effect estimate: SF - NSF
D_sen_temp = D_sen;
D_sen_temp(~tempMask) = 0;

plot_hemispheres(D_sen_temp(1,1:nHalf)',{surf_lh});
colormap(flipud(blue));
set(findall(gcf,'type','axes'),'CLim',[-1 0]);
title('Sensitivity cohort: SF vs NSF');

% 8.4 FDR significance mask: significant vertices are black
SigMap = zeros(size(D_sen));
SigMap(h_sen_temp) = -1;

plot_hemispheres(SigMap(1,1:nHalf)',{surf_lh});
colormap(flipud(gray));
set(findall(gcf,'type','axes'),'CLim',[-1 0]);
title('Sensitivity cohort: SF vs NSF, temporal PFDR < 0.05');