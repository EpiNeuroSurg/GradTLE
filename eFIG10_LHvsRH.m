clear; clc;

%% 1. Load data
load_demo
load_surface
load_color

load('/host/verges/tank/data/mo/demo/02Voxel_featsPET_hmz_flip_asym_wPDB.mat')

%% 2. Combine variables across discovery, validation, and sensitivity cohorts

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

Side_allPTs_raw = [
    Demo_disPTs.Side;
    Demo_valPTs.Side;
    Demo_senPTs.Side
];

SZFreq_allPTs = [
    Demo_disPTs.SZFreq;
    Demo_valPTs.SZFreq;
    Demo_senPTs.SZFreq
];

%% 3. Clean side and seizure-frequency variables

Side_allPTs = Side_allPTs_raw;

Side_allPTs(strcmpi(Side_allPTs_raw,'Left'))  = {'L'};
Side_allPTs(strcmpi(Side_allPTs_raw,'Right')) = {'R'};
Side_allPTs(strcmpi(Side_allPTs_raw,'left'))  = {'L'};
Side_allPTs(strcmpi(Side_allPTs_raw,'right')) = {'R'};

SZFreq_allPTs_num = zeros(size(SZFreq_allPTs));

SZFreq_allPTs_num(strcmpi(SZFreq_allPTs,'Daily'))   = 1;
SZFreq_allPTs_num(strcmpi(SZFreq_allPTs,'Weekly'))  = 2;
SZFreq_allPTs_num(strcmpi(SZFreq_allPTs,'Monthly')) = 3;
SZFreq_allPTs_num(strcmpi(SZFreq_allPTs,'Yearly'))  = 4;

%% 4. Check demographic numbers

fprintf('\n===== All patient cohorts =====\n');
fprintf('All patients: n = %d\n', numel(Side_allPTs));
fprintf('Left-sided TLE: n = %d\n', sum(strcmp(Side_allPTs,'L')));
fprintf('Right-sided TLE: n = %d\n', sum(strcmp(Side_allPTs,'R')));

%% 5. Load PET feature

MidPET_allPTs = Voxel_allPTs_neofeats_hmz_flip_asym_wPDB.MidPET;

fprintf('\nPET matrix: %d subjects x %d vertices\n', ...
    size(MidPET_allPTs,1), size(MidPET_allPTs,2));

if size(MidPET_allPTs,1) ~= numel(Side_allPTs)
    error('Subject number mismatch: PET matrix has %d rows, demographic variables have %d subjects.', ...
        size(MidPET_allPTs,1), numel(Side_allPTs));
end

%% 6. Define ipsilateral temporal cortex mask using allTempLbls

nV = size(MidPET_allPTs,2);
nHalf = nV/2;

tempMask = false(1,nV);
tempMask(1:nHalf) = ismember(economo_parc(1:nHalf), allTempLbls);

fprintf('\nTemporal labels used: %s\n', mat2str(allTempLbls));
fprintf('Temporal cortex vertices used for FDR: %d\n', sum(tempMask));

%% 7. Surface-wise GLM: left-sided vs right-sided TLE

Grp_Cov    = term(Side_allPTs);
Age_Cov    = term(Age_allPTs);
Sex_Cov    = term(Sex_allPTs);
Dur_Cov    = term(Dur_allPTs);
SZFreq_Cov = term(SZFreq_allPTs_num);

Model = 1 + Grp_Cov + Age_Cov + Sex_Cov + Dur_Cov + SZFreq_Cov;

slm_side = SurfStatLinMod(MidPET_allPTs, Model);
slm_side = SurfStatT(slm_side, Grp_Cov.L - Grp_Cov.R);

p_side = 2 * tcdf(-abs(slm_side.t), slm_side.df);
D_side = slm_side.ef;

%% 8. FDR correction within allTempLbls only

p_temp = p_side(tempMask);

[h_temp, crit_p_temp, adj_p_temp] = fdr_bh(p_temp, 0.01, 'pdep', 'yes');

h_side_temp = false(size(p_side));
h_side_temp(tempMask) = h_temp;

fprintf('\n===== Left-sided vs Right-sided TLE, allTempLbls only =====\n');
fprintf('Minimum temporal uncorrected P = %.6f\n', min(p_temp));
fprintf('Temporal FDR critical P = %.6f\n', crit_p_temp);
fprintf('Number of temporal FDR-significant vertices = %d\n', sum(h_side_temp));

%% 9. Plot maps

idxL = strcmp(Side_allPTs,'L');
idxR = strcmp(Side_allPTs,'R');

% 9.1 Left-sided TLE mean HypoPET
PET_L = mean(MidPET_allPTs(idxL,:),1);
PET_L(~tempMask) = 0;

plot_hemispheres(PET_L(1,1:nHalf)',{surf_lh});
colormap(flipud(blue));
set(findall(gcf,'type','axes'),'CLim',[-6 0]);
title('Left-sided TLE');

% 9.2 Right-sided TLE mean HypoPET
PET_R = mean(MidPET_allPTs(idxR,:),1);
PET_R(~tempMask) = 0;

plot_hemispheres(PET_R(1,1:nHalf)',{surf_lh});
colormap(flipud(blue));
set(findall(gcf,'type','axes'),'CLim',[-6 0]);
title('Right-sided TLE');

% 9.3 Raw Cohen''s d within allTempLbls
D_side_temp = D_side;
D_side_temp(~tempMask) = 0;

plot_hemispheres(D_side_temp(1,1:nHalf)',{surf_lh});
colormap(flipud(blue));
set(findall(gcf,'type','axes'),'CLim',[min(D_side_temp(1,1:nHalf)') 0]);
title('Left vs Right TLE, Cohen''s d');

% 9.4 FDR significance mask: significant vertices are black
SigMap = zeros(size(D_side));
SigMap(h_side_temp) = -1;

plot_hemispheres(SigMap(1,1:nHalf)',{surf_lh});
colormap(flipud(gray));
set(findall(gcf,'type','axes'),'CLim',[-1 0]);
title('Left vs Right TLE, temporal PFDR < 0.05');