clear;clc;

load_color
%% 1. Cohort with parliament plot
group_counts = [ 165, 62, 27, 11, 56];
group_labels = {'BTTH1', 'BFTH', 'BTTH2', 'BSBH', 'BTTH3'};
group_colors = [
                blue(192,:);
                cmap_mean(128,:);
                blue(160,:);
                cmap_orange(128,:);
                blue(128,:);
                ];

parliament_plot(group_counts, group_labels, group_colors);



%% 2. Cohort with parliament plot
group_counts = [100, 100];
group_labels = {'CERMEP','BTTH4'};
group_colors = [cmap_G2(128,:);
                
                cmap_G1(128,:)];

parliament_plot(group_counts, group_labels, group_colors);