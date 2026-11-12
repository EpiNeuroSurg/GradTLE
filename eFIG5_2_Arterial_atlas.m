clear; clc;

%% 1. Load surface and arterial atlas

load_surface

atlasDir = '/host/verges/tank/data/mo/atlas/ArterialAtlas';

fileL = fullfile( ...
    atlasDir, ...
    'lh.ArterialAtlas.32k_fsLR.func.gii');

assert(isfile(fileL), ...
    'Missing arterial atlas: %s', fileL);

gL = gifti(fileL);
labelL = round(double(gL.cdata(:)));

assert(numel(labelL) == 32492, ...
    'Expected 32492 vertices, found %d.', numel(labelL));


%% 2. Define the three arterial territories

label_ACTPL = 23;
label_MCATL = 11;
label_PCATL = 17;


%% 3. Define cortex and medial wall

[mask_lh, ~] = load_mask('midline');

cortex_mask = logical(mask_lh(:));
medial_wall = ~cortex_mask;


%% 4. Construct display data

% 0 = medial wall
% 1 = remaining cortex
% 2 = ACTPL
% 3 = MCATL
% 4 = PCATL

plotL = ones(32492,1);

plotL(medial_wall) = 0;

plotL(labelL == label_ACTPL & cortex_mask) = 2;
plotL(labelL == label_MCATL & cortex_mask) = 3;
plotL(labelL == label_PCATL & cortex_mask) = 4;

fprintf('ACTPL vertices: %d\n', nnz(plotL == 2));
fprintf('MCATL vertices: %d\n', nnz(plotL == 3));
fprintf('PCATL vertices: %d\n', nnz(plotL == 4));


%% 5. Define colors

purple_G1  = [ 93  75 183]/255;
green_skew = [102 194 165]/255;
red_mean   = [158   1  66]/255;

custom_cmap = [
    1.00 1.00 1.00       % medial wall
    0.82 0.82 0.82       % remaining cortex
                 % ACTPL
    green_skew            % MCATL
    purple_G1
    red_mean              % PCATL
    ];


%% 6. Plot using the normal plot_hemispheres appearance

plot_hemispheres( ...
    plotL, ...
    {inflated_surf_lh}, ...
    'views', 'lm');

fig = gcf;

set(fig, ...
    'Color', 'w', ...
    'Name', 'Left arterial territories', ...
    'Renderer', 'opengl');

colormap(fig, custom_cmap);

% Remove colorbar
cb = findobj(fig, 'Type', 'ColorBar');
delete(cb);

data_axes = findobj(fig, 'Type', 'axes');
set(data_axes, 'CLim', [-0.5 4.5]);


%% 7. Remove only glossy highlights

% Keep the original plot_hemispheres lighting and cortical shading.
% Do not add camlight and do not change AmbientStrength or DiffuseStrength.
brain_patch = findobj(fig, 'Type', 'patch');

set(brain_patch, ...
    'SpecularStrength', 0, ...
    'SpecularColorReflectance', 0);


%% 8. Calculate boundaries

tri = inflated_surf_lh.tri;
coord = inflated_surf_lh.coord';

edges = [
    tri(:,[1 2]);
    tri(:,[2 3]);
    tri(:,[3 1])
    ];

value1 = plotL(edges(:,1));
value2 = plotL(edges(:,2));

is_boundary = ...
    (value1 ~= value2) & ...
    ((value1 >= 2) | (value2 >= 2));

boundary_edges = unique( ...
    sort(edges(is_boundary,:),2), ...
    'rows');


%% 9. Prepare boundary coordinates

v1 = coord(boundary_edges(:,1),:);
v2 = coord(boundary_edges(:,2),:);

X = [
    v1(:,1), ...
    v2(:,1), ...
    nan(size(v1,1),1)
    ]';

Y = [
    v1(:,2), ...
    v2(:,2), ...
    nan(size(v1,1),1)
    ]';

Z = [
    v1(:,3), ...
    v2(:,3), ...
    nan(size(v1,1),1)
    ]';


%% 10. Draw black boundaries

for k = 1:numel(data_axes)

    hold(data_axes(k), 'on');

    plot3( ...
        data_axes(k), ...
        X(:), Y(:), Z(:), ...
        'Color', [0 0 0], ...
        'LineStyle', '-', ...
        'LineWidth', 0.5);

end


%% 11. Create a second figure: black whole-brain mask

% Include every surface vertex, including the medial wall
black_brain = zeros(32492,1);

plot_hemispheres( ...
    black_brain, ...
    {inflated_surf_lh}, ...
    'views', 'lm');

mask_fig = gcf;

set(mask_fig, ...
    'Color', 'w', ...
    'Name', 'Left hemisphere black mask', ...
    'Renderer', 'opengl');

% Entire surface is pure black
colormap(mask_fig, [0 0 0]);

% Remove colorbar
delete(findobj(mask_fig, 'Type', 'ColorBar'));

mask_axes = findobj(mask_fig, 'Type', 'axes');
set(mask_axes, ...
    'CLim', [-0.5 0.5], ...
    'Color', 'w');

% Remove lighting so the mask is uniformly black
delete(findall(mask_fig, 'Type', 'light'));

mask_patch = findobj(mask_fig, 'Type', 'patch');

set(mask_patch, ...
    'FaceColor', [0 0 0], ...
    'FaceLighting', 'none', ...
    'EdgeColor', 'none');

% Keep the surrounding background pure white
set(mask_fig, 'InvertHardcopy', 'off');
