clear; clc; close all;

%% ========================================================================
% eFIG8_residual_severe_HypoPET_SF_NSF.m
%
% RESIDUAL SEVERE HYPOMETABOLISM
%
% Step 1:
%   Individual severe HypoPET mask
%       MidPET w-score <= -1.96
%
% Step 2:
%   Individual residual severe HypoPET
%       SevereHypo AND NOT SurgicalCavity
%
% Step 3:
%   Compare residual severe HypoPET between:
%       SF vs NSF
%
% Statistics:
%   Vertex-wise SurfStat
%   Covariates:
%       Age
%       Sex
%       Duration
%       Seizure frequency
%
% Outputs:
%   Figure 1 = SF residual severe HypoPET
%   Figure 2 = NSF residual severe HypoPET
%   Figure 3 = SF - NSF FDR difference
%   Figure 4 = black inflated brain
%
% No saving
% No internal functions
%% ========================================================================


%% 1. LOAD
load_demo
load_surface
load_color

load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat');


%% ========================================================================
% 2. CLINICAL VARIABLES
%% ========================================================================

Age_PTs = Demo_disPTs.AgeSurg;
Sex_PTs = Demo_disPTs.Sex;
Dur_PTs = Demo_disPTs.Duration;

SZFreq_PTs = Demo_disPTs.SZFreq;

SZFreq_PTs_num = zeros(size(SZFreq_PTs));

SZFreq_PTs_num(strcmp(SZFreq_PTs,'Daily'))   = 1;
SZFreq_PTs_num(strcmp(SZFreq_PTs,'Weekly'))  = 2;
SZFreq_PTs_num(strcmp(SZFreq_PTs,'Monthly')) = 3;
SZFreq_PTs_num(strcmp(SZFreq_PTs,'Yearly'))  = 4;


%% Outcome
Prog_PTs = Demo_disPTs.x1YProgI;

Prog_bi = repmat({'NSF'},size(Prog_PTs));
Prog_bi(strcmp(Prog_PTs,'ILAE 1')) = {'SF'};

idx_SF  = strcmp(Prog_bi,'SF');
idx_NSF = strcmp(Prog_bi,'NSF');


fprintf('\n====================================================\n');
fprintf('RESIDUAL SEVERE HYPOPET ANALYSIS\n');
fprintf('====================================================\n');

fprintf('SF  = %d\n',sum(idx_SF));
fprintf('NSF = %d\n',sum(idx_NSF));


%% ========================================================================
% 3. INDIVIDUAL SEVERE HYPOPET MASK
%
% MidPET = individual PET w-score
%
% Severe HypoPET:
%       w <= -1.96
%% ========================================================================

PET_w = ...
    Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET;


nSub  = size(PET_w,1);
nVert = size(PET_w,2);
nLH   = nVert/2;


if nSub ~= height(Demo_disPTs)

    error('PET subject number does not match Demo_disPTs.');

end


%% Individual binary severe-HypoPET mask
SevereHypo = ...
    PET_w <= -1.96;


%% Remove medial wall
SevereHypo(:,economo_parc == 0) = false;


fprintf('\nThreshold = MidPET w-score <= -1.96\n');

fprintf('Mean severe vertices / patient = %.1f\n', ...
    mean(sum(SevereHypo,2)));


%% ========================================================================
% 4. LOAD INDIVIDUAL SURGICAL CAVITY
%
% Same as FIG6A_2:
%
% Right-sided patients are flipped so that
% ipsilateral hemisphere is always LEFT.
%% ========================================================================

PATH_PTs = ...
    '/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/derivatives/';


surgMask_PTs = ...
    false(nSub,nVert);


for i = 1:nSub

    File_LH = ...
        strcat( ...
        PATH_PTs, ...
        'surgMask/', ...
        Demo_disPTs.ID{i}, ...
        '/ses-post/', ...
        Demo_disPTs.ID{i}, ...
        '_ses-post_hemi-L_space-nativepro_surf-fsLR-32k_label-cavity.func.gii');


    File_RH = ...
        strcat( ...
        PATH_PTs, ...
        'surgMask/', ...
        Demo_disPTs.ID{i}, ...
        '/ses-post/', ...
        Demo_disPTs.ID{i}, ...
        '_ses-post_hemi-R_space-nativepro_surf-fsLR-32k_label-cavity.func.gii');


    if ~exist(File_LH,'file')
        error('Missing LH cavity: %s',File_LH);
    end


    if ~exist(File_RH,'file')
        error('Missing RH cavity: %s',File_RH);
    end


    gL = gifti(File_LH);
    gR = gifti(File_RH);


    surgMask_raw = ...
        [ ...
        gL.cdata', ...
        gR.cdata' ...
        ];


    %% ------------------------------------------------------------
    % Right-sided TLE:
    % swap hemispheres
    %% ------------------------------------------------------------

    if strcmp(Demo_disPTs.Side{i},'R')

        surgMask_temp = ...
            [ ...
            surgMask_raw(1,nLH+1:end), ...
            surgMask_raw(1,1:nLH) ...
            ];

    else

        surgMask_temp = ...
            surgMask_raw;

    end


    %% Binarize
    surgMask_PTs(i,:) = ...
        surgMask_temp > 0;

end


fprintf('Loaded surgical cavities = %d\n',nSub);


%% ========================================================================
% 5. INDIVIDUAL RESIDUAL SEVERE HYPOPET
%
% Residual severe HypoPET =
%
%       severe HypoPET
%       AND
%       NOT surgical cavity
%
% Each row = one patient's individual residual mask
%% ========================================================================

ResidualHypo = ...
    SevereHypo & ...
    ~surgMask_PTs;


%% Remove medial wall
ResidualHypo(:,economo_parc == 0) = false;


fprintf('\nMean residual severe vertices / patient = %.1f\n', ...
    mean(sum(ResidualHypo,2)));


%% ========================================================================
% 6. SF AND NSF RESIDUAL OCCURRENCE MAPS
%
% At each vertex:
%
% value =
% proportion of patients with residual severe HypoPET
%% ========================================================================

Residual_SF = ...
    mean( ...
        double(ResidualHypo(idx_SF,:)), ...
        1);


Residual_NSF = ...
    mean( ...
        double(ResidualHypo(idx_NSF,:)), ...
        1);


%% medial wall = 0
Residual_SF(economo_parc == 0) = 0;
Residual_NSF(economo_parc == 0) = 0;


%% ========================================================================
% 7. VERTEX-WISE SF vs NSF
%% ========================================================================

Grp_Cov    = term(Prog_bi);
Age_Cov    = term(Age_PTs);
Sex_Cov    = term(Sex_PTs);
Dur_Cov    = term(Dur_PTs);
SZFreq_Cov = term(SZFreq_PTs_num);


Model = ...
    1 + ...
    Grp_Cov + ...
    Age_Cov + ...
    Sex_Cov + ...
    Dur_Cov + ...
    SZFreq_Cov;


slm_residual = ...
    SurfStatLinMod( ...
        double(ResidualHypo), ...
        Model);


%% SF - NSF
slm_residual = ...
    SurfStatT( ...
        slm_residual, ...
        Grp_Cov.SF - Grp_Cov.NSF);


%% Same statistical direction as FIG6A_2
p_residual = ...
    tcdf( ...
        slm_residual.t, ...
        slm_residual.df);


%% Cohen's d
D_residual = ...
    slm_residual.ef;


%% ========================================================================
% 8. FDR ACROSS LABELLED CORTEX
%% ========================================================================

cortex_mask = ...
    economo_parc ~= 0;


cortex_idx = ...
    find(cortex_mask');


p_cortex = ...
    p_residual(cortex_mask');


[h_residual,crit_p_residual] = ...
    fdr_bh( ...
        p_cortex, ...
        0.05, ...
        'pdep', ...
        'yes');


sigMask_residual_FDR = ...
    false(1,nVert);


sigMask_residual_FDR( ...
    cortex_idx(h_residual == 1)) = true;


fprintf('\n----------------------------------------------------\n');
fprintf('VERTEX-WISE RESIDUAL RESULTS\n');
fprintf('----------------------------------------------------\n');

fprintf('FDR critical P       = %.8g\n', ...
    crit_p_residual);

fprintf('Significant vertices = %d\n', ...
    sum(sigMask_residual_FDR));


%% ========================================================================
% 9. PEAK FDR-SIGNIFICANT EFFECT
%% ========================================================================

if any(sigMask_residual_FDR)

    sig_vertices = ...
        find(sigMask_residual_FDR);


    D_sig = ...
        D_residual(sigMask_residual_FDR);


    [peak_abs_d,peak_local] = ...
        max(abs(D_sig));


    peak_d = ...
        D_sig(peak_local);


    peak_vertex = ...
        sig_vertices(peak_local);


    peak_p = ...
        p_residual(peak_vertex);


    fprintf('\nPeak |Cohen''s d| = %.3f\n',peak_abs_d);
    fprintf('Signed peak d     = %.3f\n',peak_d);
    fprintf('Peak vertex       = %d\n',peak_vertex);
    fprintf('Raw P at peak     = %.8g\n',peak_p);
    fprintf('P_FDR             < 0.05\n');

else

    peak_d = NaN;

    fprintf('\nNo vertices survive FDR correction.\n');

end


%% ========================================================================
% 10. COMMON COLOR RANGE FOR SF AND NSF
%% ========================================================================

Residual_max = ...
    max( ...
        [Residual_SF Residual_NSF], ...
        [], ...
        'omitnan');


Residual_max = ...
    ceil(Residual_max * 100) / 100;


if Residual_max <= 0
    Residual_max = 1;
end


fprintf('\nResidual occurrence colorbar = [0 %.2f]\n', ...
    Residual_max);


%% ========================================================================
% FIGURE 1
% SF RESIDUAL SEVERE HYPOPET
%
% Inflated ipsilateral hemisphere
%% ========================================================================

plot_hemispheres( ...
    Residual_SF(1,1:nLH)', ...
    {inflated_surf_lh});


colormap(cmap_purple);


set( ...
    findall(gcf,'type','axes'), ...
    'CLim',[0 Residual_max]);


title(sprintf( ...
    'SF residual severe HypoPET (n = %d)', ...
    sum(idx_SF)));


%% ========================================================================
% FIGURE 2
% NSF RESIDUAL SEVERE HYPOPET
%% ========================================================================

plot_hemispheres( ...
    Residual_NSF(1,1:nLH)', ...
    {inflated_surf_lh});


colormap(cmap_purple);


set( ...
    findall(gcf,'type','axes'), ...
    'CLim',[0 Residual_max]);


title(sprintf( ...
    'NSF residual severe HypoPET (n = %d)', ...
    sum(idx_NSF)));


%% ========================================================================
% 11. FDR DIFFERENCE MAP
%
% SF - NSF
%% ========================================================================

D_residual_FDR = ...
    D_residual .* ...
    sigMask_residual_FDR;


%% medial wall = 0
D_residual_FDR(economo_parc == 0) = 0;


%% Difference color range
if any(sigMask_residual_FDR)

    Dmax = ...
        max( ...
            abs(D_residual(sigMask_residual_FDR)), ...
            [], ...
            'omitnan');

else

    Dmax = ...
        max( ...
            abs(D_residual(cortex_mask')), ...
            [], ...
            'omitnan');

end


Dmax = ...
    ceil(Dmax * 100) / 100;


if Dmax <= 0
    Dmax = 1;
end


fprintf('Difference Cohen''s d colorbar = [%.2f %.2f]\n', ...
    -Dmax,Dmax);


%% ========================================================================
% FIGURE 3
% SF - NSF RESIDUAL DIFFERENCE
%% ========================================================================

plot_hemispheres( ...
    D_residual_FDR(1,1:nLH)', ...
    {inflated_surf_lh});


colormap(blue2red);


set( ...
    findall(gcf,'type','axes'), ...
    'CLim',[-Dmax Dmax]);


if any(sigMask_residual_FDR)

    title(sprintf( ...
        'Residual severe HypoPET: SF - NSF; peak d = %.3f', ...
        peak_d));

else

    title( ...
        'Residual severe HypoPET: SF - NSF; n.s.');

end


%% ========================================================================
% FIGURE 4
% WHOLE-BRAIN BLACK INFLATED SURFACE
%% ========================================================================

black_brain = ...
    zeros(nLH,1);


plot_hemispheres( ...
    black_brain, ...
    {inflated_surf_lh});


colormap([0 0 0]);


set( ...
    findall(gcf,'type','axes'), ...
    'CLim',[0 1]);


title('Whole-brain black surface');


%% ========================================================================
% 12. FINAL SUMMARY
%% ========================================================================

fprintf('\n====================================================\n');
fprintf('FINAL RESIDUAL HYPOPET SUMMARY\n');
fprintf('====================================================\n');

fprintf('Definition:\n');
fprintf('  Severe HypoPET = individual MidPET w <= -1.96\n');
fprintf('  Residual HypoPET = Severe HypoPET AND NOT cavity\n\n');

fprintf('SF  n = %d\n',sum(idx_SF));
fprintf('NSF n = %d\n',sum(idx_NSF));

fprintf('\nMean severe vertices / patient:\n');
fprintf('  SF  = %.1f\n', ...
    mean(sum(SevereHypo(idx_SF,:),2)));

fprintf('  NSF = %.1f\n', ...
    mean(sum(SevereHypo(idx_NSF,:),2)));

fprintf('\nMean residual severe vertices / patient:\n');
fprintf('  SF  = %.1f\n', ...
    mean(sum(ResidualHypo(idx_SF,:),2)));

fprintf('  NSF = %.1f\n', ...
    mean(sum(ResidualHypo(idx_NSF,:),2)));

fprintf('\nVertex-wise SF vs NSF:\n');
fprintf('  FDR-significant vertices = %d\n', ...
    sum(sigMask_residual_FDR));

if any(sigMask_residual_FDR)

    fprintf('  Peak Cohen''s d = %.3f\n', ...
        peak_d);

    fprintf('  Raw P at peak   = %.8g\n', ...
        peak_p);

    fprintf('  Peak P_FDR      < 0.05\n');

else

    fprintf('  FDR result = n.s.\n');

end

fprintf('====================================================\n');