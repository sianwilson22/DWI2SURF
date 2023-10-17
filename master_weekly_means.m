%% This script will calculate the average diffusion/surf metrics in each gestational week and 
% project them onto the surface

% Directory containing the subject surface files
main = ('/Users/sianwilson/Desktop/BCH/MASTER/');
main_input = ([main 'dHCP-fetal-final-cohort/']);
data_input = ([main_input 'dHCP-fetal-dwi-and-surf/']);
output = ([main_input 'weekly_mean_maps/']);
mkdir (output)
% Add path to functions
addpath([main 'Scripts']);
% Set number of vertices
vertices = 40962;
% Group of subjects
subject_list = importdata ([main 'Subject-Lists/age-sorted-subjects.txt']);
age_list = importdata ([main 'Subject-Lists/age-sorted-ga-rounded.txt']);
no_space_list = regexprep(subject_list, '\s', '');
% Create struct containing all surface metrics for each subject
% Define metrics being used
surf_metric = {'depth';'area';'mc'};
hemis = {'l','r'};
dwi_metric = {'tissue';'fluid';'fa';'md'};
% Define number of vertices in surface
v = 40962;
% Load in the .mat files containing all cohort values
load([data_input 'all_dwi_metric_x2.mat'])
load([data_input 'all_surf_metric.mat'])

%% Average diffusion for each gestational week
age = 25:36;
layer = {'inner';'outer'};
for hemi = 1:2
    for l = 1:2
        mean_dwi_all.(layer{l}).one.(hemis{hemi}) = cell(length(dwi_metric),1);
        mean_dwi_all.(layer{l}).two.(hemis{hemi}) = cell(length(dwi_metric),1);
        for metric = 1:4
            mean_dwi_all.(layer{l}).one.(hemis{hemi}){metric} = zeros(length(age),v);
            mean_dwi_all.(layer{l}).two.(hemis{hemi}){metric} = zeros(length(age),v);
            for g = 1:length(age)
                i = age_list() == age(g);
                % Mean dwi metrics for each vertex in each gestational week
                mean_dwi_all.(layer{l}).one.(hemis{hemi}){metric}(g,:) = mean(all_dwi_metric.(layer{l}).one.(hemis{hemi}){metric,1}(i,:));
                mean_dwi_all.(layer{l}).two.(hemis{hemi}){metric}(g,:) = mean(all_dwi_metric.(layer{l}).two.(hemis{hemi}){metric,1}(i,:));
            end
        end
    end
end
cd(data_input)
save('weekly_vertex_mean_dwi_all_metric.mat', 'mean_dwi_all');
%% Average surface metrics for each gestational week
for hemi = 1:2
    mean_surf_all.(hemis{hemi})= cell(length(surf_metric),1);
    for metric = 1:3
        mean_surf_all.(hemis{hemi}){metric} = zeros(length(age),v);
        for g = 1:length(age)
            i = age_list() == age(g);            
            % Mean surface metrics for each vertex
            mean_surf_all.(hemis{hemi}){metric}(g,:) = mean(all_surf_metric.(hemis{hemi}){metric,1}(i,:));
        end
    end
end
cd (data_input)
save('weekly_vertex_mean_surf_all_metric.mat', 'mean_surf_all');
%% DWI projected to the surface
%Example subjects for each gestational week-
surf_subjects = table2cell(readtable([main 'Subject-Lists/example_surface_subjects.csv']));

%Set the value limits for surface figure, for each dwi metric
limits = {[0.2;0.36];[0;0.1];[0.05;0.25];[0.0012;0.0018]};
colour = {(flipud(hot));(flipud(parula));(flipud(hot));(flipud(parula))};

%Stop all the figures from popping up
%set(groot,'defaultFigureVisible','on'); 

for surf = 1:length(dwi_metric)
    mkdir ([output '/' dwi_metric{surf}])
    for l = 1:2
        for age = 1 : length(ga)
            s = surf_subjects{age};
            cd ([data_input s '/surfaces/']);
            surf_rs = SurfStatReadSurf({'lh.smoothwm.native.rsl.obj',...
                'rh.smoothwm.native.rsl.obj'});
            
            cd (output);
            
            values = [mean_dwi_all.(layer{l}).two.l{surf}(age,:) ...
                mean_dwi_all.(layer{l}).two.r{surf}(age,:)];
            figure; 
            SurfStatView(values,surf_rs);
            SurfStatColLim(limits{surf});
            colormap (colour{surf});
            exportgraphics(gcf, strcat([dwi_metric{surf} '/' num2str(ga(age))...
                '_average_' dwi_metric{surf} '_map_' layer{l} '_two.jpg']),...
            'Resolution',300)
            hold off
            
        end
    end
end

%% Macrostructure means projected to the surface
surf_range = {[0, 12];[0, 1.2];[-0.6, 0.6]};
% RUN COLOURSCALE.M FIRST TO GET 'FULLSMOOTH' COLOURBAR
surf_colour ={(flipud(hot));(flipud(hot));(full_smooth)};

for surf = 2:length(surf_metric)
    mkdir ([output '/' surf_metric{surf}])
    for age = 1 : length(ga)
        s = surf_subjects{age};
        cd ([data_input s '/surfaces/']);
        surf_rs = SurfStatReadSurf({'lh.smoothwm.native.rsl.obj',...
            'rh.smoothwm.native.rsl.obj'});
        
        cd (output);
        
        values = [mean_surf_all.l{surf}(age,:) ...
            mean_surf_all.r{surf}(age,:)];
        figure;
        SurfStatView(values,surf_rs);
        SurfStatColLim(surf_range{surf});
        colormap (surf_colour{surf});
        exportgraphics(gcf, strcat([surf_metric{surf} '/' num2str(ga(age))...
            '_average_' surf_metric{surf} '_map.jpg']),...
            'Resolution',300)
        hold off
        
    end
end
