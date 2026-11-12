clear;clc;

%% 1. load data, surface, color
load_surface
load_color

%% 3. Plot cyto surface
% Create a custom colormap that maps each label to its specific color
max_label = max(allTempLbls);
custom_cmap = repmat([0.9 0.9 0.9], max_label + 1, 1);  
custom_cmap(1, :) = [1 1 1];  % White for label 0 (medial wall)
% Map labels to colors (index = label + 1 because MATLAB indexing starts at 1)

% peri -> purple
custom_cmap(12+1, :) = purple_G1;
custom_cmap(17+1, :) = purple_G1;
custom_cmap(18+1, :) = purple_G1;
% mesial -> red
custom_cmap(14+1, :) = red_mean;
custom_cmap(15+1, :) = red_mean;
custom_cmap(16+1, :) = red_mean;
% lateral -> green
custom_cmap(37+1, :) = green_skew;
custom_cmap(38+1, :) = green_skew;
custom_cmap(39+1, :) = green_skew;
custom_cmap(40+1, :) = green_skew;
custom_cmap(41+1, :) = green_skew;
custom_cmap(42+1, :) = green_skew;
custom_cmap(43+1, :) = green_skew;

obj = plot_hemispheres(economo_parc(1:end/2), {inflated_surf_lh}, 'views', 'lm');
colormap(custom_cmap);
set(findall(gcf, 'type', 'axes'), 'CLim', [0 max_label]);

%  3.1 calculate the boundary
% Extract LEFT hemisphere labels & surface
nvert = size(inflated_surf_lh.coord,2);
label_L = double(economo_parc(1:nvert));   
tri_L   = inflated_surf_lh.tri;
coord_L = inflated_surf_lh.coord';

% Build all mesh edges
edges = [tri_L(:,[1 2]);
         tri_L(:,[2 3]);
         tri_L(:,[3 1])];
lab1 = label_L(edges(:,1));
lab2 = label_L(edges(:,2));
roi1 = ismember(lab1, allTempLbls);
roi2 = ismember(lab2, allTempLbls);

is_internal = roi1 & roi2 & (lab1 ~= lab2);
is_outer = (roi1 & ~roi2) | (~roi1 & roi2);
is_boundary = is_internal | is_outer;
boundary_edges = unique(sort(edges(is_boundary,:),2), 'rows');

v1 = coord_L(boundary_edges(:,1),:);
v2 = coord_L(boundary_edges(:,2),:);

X = [v1(:,1) v2(:,1) nan(size(v1,1),1)]';
Y = [v1(:,2) v2(:,2) nan(size(v1,1),1)]';
Z = [v1(:,3) v2(:,3) nan(size(v1,1),1)]';
X = X(:); Y = Y(:); Z = Z(:);

% Draw on ALL surface axes (BrainStat makes >1 axis; exclude colorbar)
axs = findobj(gcf,'Type','axes');
axs = axs(~strcmp(get(axs,'Tag'),'Colorbar'));
for k = 1:numel(axs)
    hold(axs(k),'on');
    plot3(axs(k), X, Y, Z, 'k-', 'LineWidth', 0.1);
end