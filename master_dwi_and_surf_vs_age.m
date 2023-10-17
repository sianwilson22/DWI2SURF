%% This script will calculate the correlation between diffusion metrics and age in the CP and SP
% project them onto the surface
% Directory containing the subject surface files
main = ('/Users/sianwilson/Desktop/BCH/MASTER/');
main_input = ([main 'dHCP-fetal-final-cohort/']);
data_input = ([main_input 'dHCP-fetal-dwi-and-surf/']);
output = ([main_input 'dwi_vs_age_corr_maps/']);
surf_output = ([main_input 'surf_vs_age_corr_maps/']);
mkdir (output)
mkdir (surf_output)
% Add path to functions
addpath([main 'Scripts']);
% Set number of vertices
vertices = 40962;
% Group of subjects
subject_list = importdata ([main 'Subject-Lists/age-sorted-subjects.txt']);
age_list = importdata ([main 'Subject-Lists/age-sorted-ga.txt']);
no_space_list = regexprep(subject_list, '\s', '');
% Create struct containing all surface metrics for each subject
% Define metrics being used
surf_metric = {'depth';'area';'mc'};
hemis = {'l','r'};
dwi_metric = {'tissue';'fluid';'fa';'md'};
layer = {'inner','outer'};
% Define number of vertices in surface
v = 40962;
% Load in the .mat files containing all cohort values
load([data_input 'all_dwi_metric_x2.mat'])
load([data_input 'all_surf_metric.mat'])
% Load in the sulcal labels
l_sulci = load([main 'surface-template-28w-labels/lh.28w.sulcal']);
r_sulci = load([main 'surface-template-28w-labels/rh.28w.sulcal']);
both_sulci = [l_sulci; r_sulci];
%% Correlation of diffusion metrics with age, projected onto the surface
%Read in surface files
cd ([data_input 'sub-CC00974XX18_ses-34631/surfaces/']); 
surf_rs = SurfStatReadSurf({'lh.smoothwm.native.rsl.obj','rh.smoothwm.native.rsl.obj'});

%Initialise matrix for diffusion correlation values
for n = 1:length(dwi_metric)
    for l = 1:2
        dwi_age_corr.(layer{l}).l = cell (length(dwi_metric),1);
        dwi_age_corr.(layer{l}).r = cell (length(dwi_metric),1); 
        for hemi = 1:2
            dwi_age_corr.(layer{l}).(hemis{hemi}){n} = zeros(length(l_sulci),1);
            for vertex = 1:length(l_sulci)
                dwi_age_corr.(layer{l}).(hemis{hemi}){n}(vertex) = corr...
                    (all_dwi_metric.(layer{l}).two.(hemis{hemi}){n}(:,vertex),age_list);
            end
        end
        dwi_age_corr_both.(layer{l}){n} = [dwi_age_corr.(layer{l}).l{n};dwi_age_corr.(layer{l}).r{n}];
        cd (output)
        figure; SurfStatView(dwi_age_corr_both.(layer{l}){n},surf_rs);
        SurfStatColLim([-0.8, 0.8])
        colormap(full_smooth);
        exportgraphics(gcf, strcat([dwi_metric{n} '_' layer{l} '_two_age_correlation.jpg']),'Resolution',300)
        hold off
    end
end

%% Correlate surface metrics with age and project onto the surface 
%Read in surface files
cd ([data_input 'sub-CC00974XX18_ses-34631/surfaces/']); 
surf_rs = SurfStatReadSurf({'lh.smoothwm.native.rsl.obj','rh.smoothwm.native.rsl.obj'});

%Initialise matrix for diffusion correlation values

for n = 1:length(surf_metric)
    for hemi = 1:2
        surf_age_corr.(hemis{hemi}) = cell(length(surf_metric),1);
    
    % Find correlation between surface metrics and age for each vertex
    for vertex = 1:length(l_sulci)
        surf_age_corr.(hemis{hemi}){n}(vertex) = corr(all_surf_metric.(hemis{hemi}){n}(:,vertex),age_list);
    end
    end
    
    %Project correlation coefficients onto the surface to see regional
    %variations
    cd (surf_output)
    figure; SurfStatView([surf_age_corr.l{n} surf_age_corr.r{n}],surf_rs);
    SurfStatColLim([0, 1])
    colormap(full_smooth)
    exportgraphics(gcf, strcat([surf_metric{n} '_age_correlation.jpg']),'Resolution',300)
    hold off
end
%% Graphs of average values across whole brain against gestational age of each subject
j  = (23:1:36);
whole_brain_mean.sp.l = cell(length(dwi_metric),1);
whole_brain_mean.cp.lh = cell(length(dwi_metric),1);

load ('/Users/sianwilson/Desktop/BCH/Surface-diffusion-study/Scripts/HJ_Smap.mat');

for d = 1:length(dwi_metric)
    for hemi = 1:2
        whole_brain_mean.sp.(hemis{hemi}) = cell(length(dwi_metric),1);
        whole_brain_mean.cp.(hemis{hemi}) = cell(length(dwi_metric),1);
        
        for subject = 1:length(age_list)
            sp_dwi.(hemis{hemi}) = [all_dwi_metric.inner.one.(hemis{hemi}){d}(subject,:)...
                all_dwi_metric.inner.two.(hemis{hemi}){d}(subject,:)];
            cp_dwi.(hemis{hemi}) = [all_dwi_metric.outer.one.(hemis{hemi}){d}(subject,:)...
                all_dwi_metric.outer.two.(hemis{hemi}){d}(subject,:)];
            
            whole_brain_mean.sp.(hemis{hemi}){d,1}(subject,:) = mean(sp_dwi.(hemis{hemi}));
            whole_brain_mean.cp.(hemis{hemi}){d,1}(subject,:) = mean(cp_dwi.(hemis{hemi}));
            
        end
    end
    
    % Evaluate the fitted polynomial p and plot:
    x1 = [age_list age_list];
    y1 = [whole_brain_mean.cp.l{d,1}(:) whole_brain_mean.cp.r{d,1}(:)];
    
    x2 = [age_list age_list];
    y2 = [whole_brain_mean.sp.l{d,1}(:) whole_brain_mean.sp.r{d,1}(:)];
    
    % Polynomial fitting
    p = polyfit(x1,y1,2);
    f1 = polyval(p,x1);
    
    p = polyfit(x2,y2,2);
    f2 = polyval(p,x2);
    
    plot_title_cp = 'Cortical plate';
    plot_title_sp = 'Subplate';
    c = Smap(d+1,:);
    
    figure;
    set(gcf,'Color','w');
    scatter(age_list,whole_brain_mean.cp.l{d,1}(:),'o','MarkerEdgeColor',...
        "#000000",'MarkerFaceColor',"#000000",'LineWidth',1.5,'MarkerFaceAlpha',0.3)
    hold on
    scatter(age_list,whole_brain_mean.cp.r{d,1}(:),'o','MarkerEdgeColor',...
        c,'MarkerFaceColor',c,'LineWidth',1.5,'MarkerFaceAlpha',0.6)
    hold on
    plot(x1,f1,'k:','LineWidth',2.5)
    title(plot_title_cp,'FontSize',16)
    xlabel('Gestational Age','FontSize',16)
    xlim([22 38]);
    ylabel ([dwi_metric{d}],'FontSize',16)
    ax = gca;
    ax.FontSize = 14;
    hold off
    exportgraphics(gcf, strcat([output '/' dwi_metric{d} '_metric_vs_age_cp.jpg']),'Resolution',300)
    
    figure;
    set(gcf,'Color','w');
    scatter(age_list,whole_brain_mean.sp.l{d,1}(:),'o','MarkerEdgeColor',...
        "#000000",'MarkerFaceColor',"#000000",'LineWidth',1.5,'MarkerFaceAlpha',0.3)
    hold on
    scatter(age_list,whole_brain_mean.sp.r{d,1}(:),'o','MarkerEdgeColor',...
        c,'MarkerFaceColor',c,'LineWidth',1.5,'MarkerFaceAlpha',0.6)
    hold on
    plot(x2,f2,'k:','LineWidth',2.5)
    title(plot_title_sp,'FontSize',16)
    xlabel('Gestational Age','FontSize',16)
    xlim([22 38]);
    ylabel ([dwi_metric{d}],'FontSize',16)
    ax = gca;
    ax.FontSize = 14;
    hold off
    exportgraphics(gcf, strcat([output '/' dwi_metric{d} '_metric_vs_age_sp.jpg']),'Resolution',300)
end

