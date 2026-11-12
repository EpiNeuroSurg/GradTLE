disp('Color definition')
blue = cbrewer('seq', 'Blues', 256); 
red = (cbrewer('seq', 'Reds', 256)); 
gray = cbrewer('seq', 'Greys', 256);
blue0black = blue; blue0black(1,:) = [0 0 0]; 
white2black = [1 1 1; 0 0 0];
blue2red = flipud([flipud(blue);red]);

colors_cyto = [...    
    126 39  126;
     48 105 156;
    167 211 140;
    254 205 7;
    243 236 25];
colors_cyto_nor = colors_cyto/255;

purple_G1   = [ 93  75 183]/255;
Gray_G2     = [120 130 150]/255;
red_mean    = [158   1  66]/255;
green_skew  = [102 194 165]/255;

% ---- Build sequential colormaps: white -> color ----
n = 256;
t = linspace(0,1,n)';% 0=low, 1=high
white = ones(n,3);

cmap_G1     = white.*(1-t) + t.*purple_G1;
cmap_G2     = white.*(1-t) + t.*Gray_G2;
cmap_mean   = white.*(1-t) + t.*red_mean;
cmap_skew   = white.*(1-t) + t.*green_skew;

% bulid for hippo
cmap_purple = white.*(1-t) + t.*(colors_cyto_nor(1,:));
cmap_blue   = white.*(1-t) + t.*(colors_cyto_nor(2,:));
cmap_green  = white.*(1-t) + t.*(colors_cyto_nor(3,:));
cmap_orange = white.*(1-t) + t.*(colors_cyto_nor(4,:));
cmap_yellow = white.*(1-t) + t.*(colors_cyto_nor(5,:));
cmap_red    = white.*(1-t) + t.*(red);

cmap_gray2red = [gray(16,:)].*(1-t) + t.*red_mean; 