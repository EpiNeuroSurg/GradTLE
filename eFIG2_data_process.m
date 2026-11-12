% clear;clc;

%% 1. Load data, surface, color
load_surface
load_color

load('/host/verges/tank/data/mo/demo/01Voxel_data_hmz.mat')
load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat')

% indi raw PET   
% run Step0_1*.m
indiRawPET = Raw_PTs_MidPET(1,:)';
indiRawPET(economo_parc == 0) = max(indiRawPET);
plot_hemispheres(indiRawPET, {surf_lh, surf_rh}, 'views', 'lm'); colormap(flipud(blue))
 
% indi hmz PET
indiHmzPET = Voxel_all_MidPET_hmz(38,:)';
indiHmzPET(economo_parc == 0) = max(indiHmzPET);
plot_hemispheres(indiHmzPET, {surf_lh, surf_rh}, 'views', 'lm'); colormap(flipud(blue))

% indi fliped PET
indiFlipPET = Voxel_PTs_neofeats_hmz_flip(1,:)';
indiFlipPET(economo_parc == 0) = max(indiFlipPET);
plot_hemispheres(indiFlipPET, {surf_lh, surf_rh}, 'views', 'lm'); colormap(flipud(blue))

% indi w-score PET
indiWscorePET = Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET(1,:)';
indiWscorePET(economo_parc == 0) = max(indiWscorePET);
plot_hemispheres(indiWscorePET, {surf_lh, surf_rh}, 'views', 'lm'); colormap(flipud(blue))