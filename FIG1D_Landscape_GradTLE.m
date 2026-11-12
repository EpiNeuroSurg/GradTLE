clear;clc;

%% 1. Load data, surface, color
load_demo
load_surface
load_color

%% 2. plot GradTLE
load('/host/verges/tank/data/mo/demo/02Parcel_feats_hmz_flip_asym_wPDB.mat')
load('/host/verges/tank/data/mo/demo/02Parcel2_feats_hmz_flip_asym_wPDB.mat')
load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat')

%  2.1 parcel-wise; Economo
Parcel_MidPET_hmz_flip_asym_wPDB_mean = mean(Parcel_PTs_feats_hmz_flip_asym_wPDB.MidPET, 1);
economo_parc_lh = economo_parc(1:end/2, 1);
vertex_data = zeros(length(economo_parc_lh), 1);
for i = 1:43
    vertex_data(economo_parc_lh == i) = Parcel_MidPET_hmz_flip_asym_wPDB_mean(i);
end
plot_hemispheres(vertex_data, surf_lh, 'views', 'lm');
colormap(flipud(blue));
set(findall(gcf,'type','axes'),'CLim',[-4 0]);

%  2.2 parcel-wise; Julich
Parcel2_MidPET_hmz_flip_asym_wPDB_mean = mean(Parcel2_PTs_feats_hmz_flip_asym_wPDB.MidPET, 1);
vertex_data = zeros(length(julich_parc(1:end/2, 1)), 1);
for i = 1:148
    vertex_data(julich_parc1000 == i) = Parcel2_MidPET_hmz_flip_asym_wPDB_mean(i);
end
vertex_data(ismember(julich_parc1000,[148 147 137 70 69 18])) = max(Parcel2_MidPET_hmz_flip_asym_wPDB_mean);
plot_hemispheres(vertex_data, surf_lh, 'views', 'lm');
colormap(flipud(blue));
set(findall(gcf,'type','axes'),'CLim',[-4 0]);

%  2.3 voxel-wise
Voxel_MidPET_hmz_flip_asym_wPDB_mean = mean(Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET, 1);
Voxel_MidPET_hmz_flip_asym_wPDB_mean(find(economo_parc == 0)) = 0;
plot_hemispheres(Voxel_MidPET_hmz_flip_asym_wPDB_mean(1,1:end/2)', surf_lh, 'views', 'lm'); 
colormap(flipud(blue))
set(findall(gcf,'type','axes'),'CLim',[-4 0]);