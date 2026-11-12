clear;clc;

%% 1. Load data, surface, color
load_demo
load_surface
load_color

%% 2. Analysis surface
%  2.3 Blank brain for distance
blank_brain = zeros(32492,1);
plot_hemispheres(blank_brain, surf_lh, 'views', 'l');
colormap([0.9 0.9 0.9])

%  2.3 Blank brain for distance
blank_brain = zeros(32492,1);
plot_hemispheres(blank_brain, surf_rh, 'views', 'l');
colormap([0 0 0])

%  2.4 Economo-atlas
economo_parc_5layers(economo_parc_5layers == 0) = NaN;
plot_hemispheres(economo_parc_5layers(1:32492,1), surf_lh, 'views', 'l');
colormap(colors_cyto_nor)

%  2.5 Background Hippo and data
bkg_hipp = zeros(16,1);
figure;
plot_subcortical(bkg_hipp, 'color_range', [0 1], 'cmap', 'gray');

load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat')
hippo_pet_419 = mean(Voxel_PTs_subfeats_hmz_flip_asym_wPDB.HippoPET(:,1:419))';
hippo_data_3784_lh = interp1(1:419, hippo_pet_419, linspace(1, 419, 3784))';
data_16 = zeros(16,1); data_16(4) = 1;
figure;
plot_subcortical_modified(data_16, 'color_range', [-4 4], 'cmap', 'RdBu_r',  'hippo_lh_data', hippo_data_3784_lh)

%  2.6 surgical cavity
assume_mask = (ismember(economo_parc, mesTempLbls) | ...
               ismember(economo_parc, latTempLbls)) & ...
              ~ismember(economo_parc, [37 41]);
assume_mask_out = ~assume_mask;

assume_surg = economo_parc(1:end/2, 1);
assume_surg(assume_mask_out(1:end/2, 1)) = 0.1;
assume_surg(assume_mask(1:end/2, 1)) = 1;
plot_hemispheres(assume_surg, surf_lh, 'views', 'l'); colormap(white2black);
% colorbar off;

