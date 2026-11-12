
close all;clear;clc;

%% 1. Load data, surface, color
load_demo
load_surface      % must provide surf_rh (same as surf_lh, for the right hemisphere)
load_color

load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat')

%% 2. Right hemisphere, vertex-wise results
Voxel_MidPET_hmz_flip_asym_wPDB_mean = mean(Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET, 1);
Voxel_MidPET_hmz_flip_asym_wPDB_mean(find(economo_parc == 0)) = 0;   % mask medial wall (both hemispheres)

n_half = length(Voxel_MidPET_hmz_flip_asym_wPDB_mean) / 2;
% The right half of the vector is empty (all zeros), so copy the LEFT hemisphere
% values onto the right hemisphere (vertex i in lh <-> vertex i in rh).
vertex_data_rh = Voxel_MidPET_hmz_flip_asym_wPDB_mean(1, 1:n_half)';       % first half = left hemisphere

plot_hemispheres(vertex_data_rh, surf_rh, 'views', 'ml');   % lateral first, then medial
colormap(flipud(blue));
set(findall(gcf,'type','axes'),'CLim',[-6 1]);

%% 3. Right hemisphere, no data (blank surface, white medial wall)
medial_wall = (economo_parc(1:end/2, 1) == 0);           % medial wall mask (lh mask, vertex-matched to rh)

vertex_data_blank = ones(size(surf_rh.coord, 2), 1);     % 1 = cortex (gray)
vertex_data_blank(medial_wall) = 0;                      % 0 = medial wall (white)

plot_hemispheres(vertex_data_blank, surf_rh, 'views', 'ml');
colormap([1 1 1; 0.9 0.9 0.9]);                          % row 1 = white (medial wall), row 2 = gray (cortex)
set(findall(gcf,'type','axes'),'CLim',[0 1]);            % colorbar now shows exactly these two colors

%% 4. Right hemisphere, whole brain black
vertex_data_black = zeros(size(surf_rh.coord, 2), 1);
plot_hemispheres(vertex_data_black, surf_rh, 'views', 'ml');
colormap([0 0 0]);

% Notes:
% - If your surf_rh is stored as surf_rh.coord (3 x N), the size() call above works.
%   If it is N x 3 (surf_rh.vertices), use: zeros(size(surf_rh.vertices,1),1)
% - If lateral/medial appear swapped for the right hemisphere, use 'ml' instead of 'lm'.