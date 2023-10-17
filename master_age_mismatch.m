%% Age mismatched coupling - Older Sulcal Depth Mean vs. younger subjects tissue fraction
% RUN COLOURSCALE.M FIRST

% Variables to change:
% Which gestational week the sulcal depth is high enough (depth average >= 3)
% 28w = 6, 31w = 9, 35w = 13
week = 13;
% Index of the label of interest (i.e Central Sulcus = 2, Superior temporal
% sulcus = 9, Inferior frontal sulcus = 5, Post-central sulcus = 7)
sulcal_label = 9;
% Which example subject to use to project values (see variable
% surf_subjects)
eg_subj = 2;
% Which dwi metric (1 = tissue, 2 = fluid, 3 = fa, 4 = md)
dwi = 1;
% Which surface metric (1 = depth, 2 = area, 3 = curvature)
sm = 1;
% Define layer subplate (inner = 1) or cortical plate (outer = 2)
l = 2;
% Define number of subjects to look at (67 = all subjects up to age 30 GW)
num = 67;

%Stop all the figures from popping up
set(groot,'defaultFigureVisible','on'); 
load('RWB.mat')
% Fixed set up variables 
main = ('/Users/sianwilson/Desktop/BCH/MASTER/');
main_input = ([main 'dHCP-fetal-final-cohort/']);
data_input = ([main_input 'dHCP-fetal-dwi-and-surf/']);
output = ([main_input 'age-mismatch-mean/']);
mkdir (output)
coupling_dir = ([main_input 'coupling/']);
% Changeable variables: define metrics being used
surf_metric = {'depth';'area';'mc'};
dwi_metric = {'tissue';'fluid';'fa';'md'};
% Degree of neighbours 
k=5;
% Define number of vertices in surface
v = 40962;
load([coupling_dir 'cat_format_lh_rh.' num2str(k) 'neighbors_of_neighbours.mat']);
% Define the hemisphere names & layers
hemis = {'l', 'r'};
layer = {'inner';'outer'};
% Load variable 'mean_surf_all_l' for mean in each vertex
load([data_input 'weekly_vertex_mean_surf_all_metric.mat']);
% Load variable 'all_dwi_metric' for individual subject dwi values
load([data_input 'all_dwi_metric_x2.mat'])
load([main 'Subject-Lists/example_biweekly_surf_subjects.mat']);
% Load subjects
subject_list = importdata ([main 'Subject-Lists/age-sorted-subjects.txt']);
age_list = importdata ([main 'Subject-Lists/age-sorted-ga-rounded.txt']);
age_biweekly = importdata ([main 'Subject-Lists/age-sorted-biweekly-groups.txt']);

%  Sulcus indexing
load([main 'surface-template-28w-labels/lh.28w.sulcal']);
load([main 'surface-template-28w-labels/rh.28w.sulcal']);
both_sulci = [lh_28w; rh_28w];
sulci = table2cell(readtable([main 'surface-template-28w-labels/sulci_label_idx.csv']));

sulcus_index.l = cell(19,1);
sulcus_index.r = cell(19,1);

for label = 1:range(both_sulci)
    % Initialise cells for the indices of major sulci and their sulcal depth
    sulcus_index.l{label} = find(lh_28w==label);
    sulcus_index.r{label} = find(rh_28w==label);
end
%% Mismatch analysis 

% Initialise structs for the mean values 
surf_mean.l = mean_surf_all.l{dwi}(week,:);
surf_mean.r = mean_surf_all.r{dwi}(week,:);
all_corr_values.l = zeros(num,length(neighbours_cat.l));
all_p_values.l = zeros(num,length(neighbours_cat.r));
all_corr_values.r = zeros(num,length(neighbours_cat.r));
all_p_values.r = zeros(num,length(neighbours_cat.r));

% Calculate mean surf values in the week of interest    
% Choose which subject to use as example
for subj_no = 1:num
    for hemi = 1:2
        cs_neighvec.([hemis{hemi}]) = neighbours_cat.([hemis{hemi}])(sulcus_index.([hemis{hemi}]){sulcal_label},:);
        subj_dwi = all_dwi_metric.([layer{l}]).two.([hemis{hemi}]){dwi}(subj_no,:);
        
        result_corr.([hemis{hemi}]) = zeros(length(cs_neighvec.([hemis{hemi}])),1);
        result_p.([hemis{hemi}]) = zeros(length(cs_neighvec.([hemis{hemi}])),1);
        
        for row = 1: length(cs_neighvec.([hemis{hemi}]))
            neighvec_input = cat(1, sulcus_index.([hemis{hemi}]){sulcal_label}(row), ...
                cs_neighvec.([hemis{hemi}]){row,:});
            m1 = surf_mean.([hemis{hemi}])(neighvec_input);
            m2 = subj_dwi(neighvec_input);
            [correlation,p] = corrcoef(m1, m2);
            result_corr.([hemis{hemi}])(row) = correlation(1,2);
            result_p.([hemis{hemi}])(row) = p(1,2);            
        end
    end  
    % Make matrix of results
    brain.l = zeros(1,v);
    brain.r = zeros(1,v);
    brain_p.l = zeros(1,v);
    brain_p.r = zeros(1,v);
    
    for hemi = 1:2      
        % Find the relevant values of sulcal depth for each subject from the
        % large struct of all depth values (left and right)
        for i = 1:length(sulcus_index.([hemis{hemi}]){sulcal_label})
            label = sulcus_index.([hemis{hemi}]){sulcal_label}(i);
            brain.([hemis{hemi}])(1,label) = result_corr.([hemis{hemi}])(i);
            brain_p.([hemis{hemi}])(1,label) = result_p.([hemis{hemi}])(i);            
        end
        all_corr_values.([hemis{hemi}])(subj_no,:) = brain.([hemis{hemi}]);
        all_p_values.([hemis{hemi}])(subj_no,:) = brain_p.([hemis{hemi}]);
    end    
    % Visualise result
    % Note - run scropt colourscale.m for 'fullsmooth' colour bar to work    
%     s = surf_subjects{eg_subj};
%     cd ([data_input s '/surfaces/']);
%     surf_rs = SurfStatReadSurf({'lh.smoothwm.native.rsl.obj','rh.smoothwm.native.rsl.obj'});
%     cd (output);
%     depth = [brain.l brain.r];
%     figure; SurfStatView(depth,surf_rs);
%     SurfStatColLim([-1 1])
%     %SurfStatColLim([1 8])
%     colormap (flipud(hot))
%     colormap(full_smooth)
%     exportgraphics(gcf, strcat([sulci{sulcal_label,2} '-' layer{l} '-' num2str(age_list(subj_no)) '-' cell2mat(subject_list(subj_no)) '_28w_mistmatch_tissue_depth_corr.jpg']),'Resolution',300)
%     hold off
end 

%% Average r values for biweekly age groups
age = 25:30;
q = 0.05;
for hemi = 1:2
    % Initialise    
    mean_pearsonr_values.([hemis{hemi}]) = zeros (length(age),length(neighbours_cat.([hemis{hemi}]))); 
    sig_mean_pearsonr_values.([hemis{hemi}]) = zeros (length(age),length(neighbours_cat.([hemis{hemi}]))); 
    test.([hemis{hemi}]).p = zeros(length(age),length(neighbours_cat.([hemis{hemi}])));
    test.([hemis{hemi}]).h = zeros(length(age),length(neighbours_cat.([hemis{hemi}])));
 
    % Run through each gestational week and calculate mean Pearson r
    for g = 1:length(age)
    i = age_list(1:num) == age(g);
    mean_pearsonr_values.([hemis{hemi}])(g,:) = trimmean(all_corr_values.([hemis{hemi}])(i,:),40);
    
%     % Single sample ttest for the correlation values in each week
%     [h,p] = ttest(all_corr_values.([hemis{hemi}])(i,:));
%     test.([hemis{hemi}]).p(g,:) = p;
%     test.([hemis{hemi}]).h(g,:) = h;
%     
%     % FDR correction - to establish which weeks have significant
%     % correlations (will be indexed in sig_weeks)
%     [pid,pn] = FDR(test.([hemis{hemi}]).p,q);
%     pID.([hemis{hemi}])(g) =  pid;
%     pN.([hemis{hemi}])(g) = pn;
%     
%     % Find out which vertices did not have p value below adjusted FDR threshold
%     idx_not_sig = (test.([hemis{hemi}]).p(g,:) >= pid);
%     
%     % Make those vertices equal to zero
%     w = mean_pearsonr_values.([hemis{hemi}])(g,:);
%     w(idx_not_sig) = 0;
%     
%     % Make new struct, just containing significant correlation values, and
%     % not significant values are zeroed
%     sig_mean_pearsonr_values.([hemis{hemi}])(g,:) = w;
%     
    end
end

%% Project pearson r values to the surface pre and post fdr correction
mkdir ([output 'images']);
for g = 1:length(age)
    s = surf_subjects{eg_subj};
    cd ([data_input s '/surfaces/']);
    surf_rs = SurfStatReadSurf({'lh.smoothwm.native.rsl.obj','rh.smoothwm.native.rsl.obj'});
   %Pre-FDR correction
    pearsonr = [mean_pearsonr_values.l(g,:) mean_pearsonr_values.r(g,:)];
    figure; SurfStatView(pearsonr,surf_rs);
    colormap(RWB)
    SurfStatColLim([-1, 1])
    exportgraphics(gcf, strcat([output 'images/' sulci{sulcal_label,2} '_' num2str(age(g))...
        '_' layer{l} '2' dwi_metric{dwi} '_vs_' surf_metric{sm} '_average_corrcoef_prefdr_map_weekly_med_RWB.jpg']),'Resolution',300)
    hold off
   %Post-FDR correction
%     sigpearsonr = [sig_mean_pearsonr_values.l(g,:) sig_mean_pearsonr_values.r(g,:)];
%     figure; SurfStatView(sigpearsonr,surf_rs);
%     colormap(full_smooth)
%     SurfStatColLim([-1, 1])
%     exportgraphics(gcf, strcat([output 'images/' sulci{sulcal_label,2} '_' num2str(age(g))...
%         '_' layer{l} '2' dwi_metric{dwi} '_vs_' surf_metric{sm} '_average_corrcoef_postfdr_map_weekly_med.jpg']),'Resolution',300)
%     hold off
end


% %% Age histogram to show what range is being analysed
% figure;
% set(gca, 'Color', 'None')
% histogram (age_list(1:80),'BinWidth',1,'FaceColor',	"#000000",'FaceAlpha','0.4')
% %%
% load([main 'Analysis_Input/RWB.mat']);
% 
% % Or look manually
% l_max_idx = 453;
% 
% neighvec_input = cat(1, sulcus_index.r{sulcal_label}(l_max_idx), cs_neighvec.r{l_max_idx,:});
% cs_m1 = surf_mean.r(neighvec_input);
% cs_m2 = subj_dwi(neighvec_input);
% 
% % Test one neighbourhood
% figure; 
% set(gcf,'Color','w');
% ax = gca;
% ax.FontSize = 16;
% hold on
% plot(cs_m1, cs_m2,'o','LineStyle', 'none',...
%     'marker','o','MarkerSize',4,'MarkerEdgeColor',RWB(4,:),'MarkerFaceColor',RWB(8,:),'LineWidth',5)
% hold on
% p = polyfit(cs_m1,cs_m2,1);
% f = polyval(p,cs_m1);
% best_fit = plot(cs_m1,f,'-r','LineWidth',1);
% 
% corrcoef(cs_m1,cs_m2)
