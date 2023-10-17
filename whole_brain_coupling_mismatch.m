%% NEIGHBOURHOOD COUPLING BETWEEN SURFACE AND DIFFUSION METRICS
% Within-subject analysis of the correlation between surface and diffusion
% values in patches of vertices

% Changeable variables: define metrics being used
surf_metric = {'depth';'area';'mc'};
dwi_metric = {'tissue';'fluid';'fa';'md'};
layer = {'inner';'outer'};
% Index for which SURF metric (1 = depth, 2 = area, 3 = curvature)
sm = 1;
% Index for which DWI metric (1 = tissue, 2 = fluid, 3 = fa, 4 = md)
dwi = 1; 
% Define layer subplate (inner = 1) or cortical plate (outer = 2)
l = 2;

% Directories containing the subject surface files, 
%neighbourhoods and outputs
main = ('/Users/sianwilson/Desktop/BCH/MASTER/');
main_input = ([main 'dHCP-fetal-final-cohort/']);
data_input = ([main_input 'dHCP-fetal-dwi-and-surf/']);
coupling_dir = ([main_input 'coupling/']);
mkdir (coupling_dir);
output = ([main_input 'whole-brain-coupling-mismatch/']);
mkdir (output);
% Add path to functions
addpath([main 'Scripts']);
% Set number of vertices
vertices = 40962;
% Group of subjects
subject_list = importdata ([main 'Subject-Lists/age-sorted-subjects.txt']);
age_list = importdata ([main 'Subject-Lists/age-sorted-ga-rounded.txt']);
age_biweekly = importdata([main 'Subject-Lists/age-sorted-biweekly-groups.txt']);
load([main 'Subject-Lists/example_biweekly_surf_subjects.mat']);
no_space_list = regexprep(subject_list, '\s', '');
% Define number of vertices in surface
v = 40962;
% Load in the .mat files containing all cohort values
load([data_input 'all_dwi_metric_x2.mat'])
load([data_input 'all_surf_metric.mat'])

% Define parameters used for neighbourhood (check master_neighbourhood)
k = 5;
fwhm = 10;
%  Load in neighbour information- variable called 'neighbours_cat'
load([coupling_dir 'cat_format_lh_rh.' num2str(k) 'neighbors_of_neighbours.mat']);
load([main 'Scripts/RWB.mat']);

%% Extra bit to load mean surface values
% Load variable 'mean_surf_all_l' for mean in each vertex
load([data_input 'weekly_vertex_mean_surf_all_metric.mat']);
week = 12;
dwi = 1;
surf_mean.l = mean_surf_all.l{dwi}(week,:);
surf_mean.r = mean_surf_all.r{dwi}(week,:);

%% Neighbourhood correlations
hemis = {'l','r'};

for hemi = 1:2
    fprintf('%s\n',hemis{hemi})
    all_corr_values.([hemis{hemi}]) = zeros(length(subject_list),...
        length(neighbours_cat.([hemis{hemi}])));
    all_p_values.([hemis{hemi}]) = zeros(length(subject_list),...
        length(neighbours_cat.([hemis{hemi}])));
    
    % Reformat neighbours list
    neigh_cat = arrayfun(@(row) vertcat(neighbours_cat.([hemis{hemi}]){row,:}), ...
        1:size(neighbours_cat.([hemis{hemi}]),1), 'UniformOutput', false);
        
    for s = 1:length(no_space_list)
        
        % Define the subject
        x = string(no_space_list{s});
        subject = regexprep(x, '\s', '' );
        % Print out name of subject to terminal
        fprintf('%s\n',subject)
        fprintf('%s\n', 'Correlation between micro and macro...');
        
        % CORRCOEF UNWEIGHTED (Pearson r)
        result_corr = zeros(length(neigh_cat),1);
        result_p = zeros(length(neigh_cat),1);
        % x value
        variables.m1 = surf_mean.([hemis{hemi}]);
        % y value
        variables.m2 = all_dwi_metric.([layer{l}]).two.([hemis{hemi}]){dwi}(s,:);
        
        for row = 1: length(neighbours_cat.l)
            neighvec_input = [row ; neigh_cat{row}];
            m1 = variables.m1(neighvec_input);
            m2 = variables.m2(neighvec_input);
            [correlation,p] = corrcoef(m1, m2);
            result_corr(row) = correlation(1,2);
            result_p(row) = p(1,2);
        end
        
        all_corr_values.([hemis{hemi}])(s,:) = result_corr;
        all_p_values.([hemis{hemi}])(s,:) = result_p;
        
        save(fullfile(output, [char(subject),'-',hemis{hemi},...
            '.k-',num2str(k), '-fwhm-',num2str(fwhm) layer{l} '2' ...
            dwi_metric{dwi} surf_metric{sm} '-corrcoef.txt']), 'result_corr','-ascii','-double');
        save(fullfile(output, [char(subject),'-',hemis{hemi}, ...
            '.k-',num2str(k), '-fwhm-',num2str(fwhm) layer{l} '2' ...
            dwi_metric{dwi} surf_metric{sm} '-corrcoef-p.txt']), 'result_p','-ascii','-double');
    end
end
corr_filename = (['k-',num2str(k), '-', num2str(24+week),'-fwhm-',num2str(fwhm) '.' layer{l} ...
    dwi_metric{dwi} '-' surf_metric{sm} '-corrcoef.mat']);
save(fullfile(output, corr_filename), 'all_corr_values','-v7');

p_filename = (['k-',num2str(k), '-', num2str(24+week),'-fwhm-',num2str(fwhm) '.' layer{l} ...
    dwi_metric{dwi} '-' surf_metric{sm} '-p.mat']);
save(fullfile(output, p_filename), 'all_p_values','-v7');

%% Weekly average Pearson r values at each vertex, single sample ttest in each gestational week for significance
age = 25:2:35;
q = 0.05;
for hemi = 1:2
    % Initialise    
    mean_pearsonr_values.([hemis{hemi}]) = zeros (length(age),length(neigh_cat)); 
    sig_mean_pearsonr_values.([hemis{hemi}]) = zeros (length(age),length(neigh_cat)); 
    test.([hemis{hemi}]).p = zeros(length(age),length(neigh_cat));
    test.([hemis{hemi}]).h = zeros(length(age),length(neigh_cat));
 
    % Run through each gestational week and calculate mean Pearson r
    for g = 1:length(age)
    i = age_biweekly() == age(g);
    mean_pearsonr_values.([hemis{hemi}])(g,:) = mean(all_corr_values.([hemis{hemi}])(i,:));
    
    % Single sample ttest for the correlation values in each week
    [h,p] = ttest(all_corr_values.([hemis{hemi}])(i,:));
    test.([hemis{hemi}]).p(g,:) = p;
    test.([hemis{hemi}]).h(g,:) = h;
    
    % FDR correction - to establish which weeks have significant
    % correlations (will be indexed in sig_weeks)
    [pid,pn] = FDR(test.([hemis{hemi}]).p,q);
    pID.([hemis{hemi}])(g) =  pid;
    pN.([hemis{hemi}])(g) = pn;
    
    % Find out which vertices did not have p value below adjusted FDR threshold
    idx_not_sig = (test.([hemis{hemi}]).p(g,:) >= pid);
    
    % Make those vertices equal to zero
    w = mean_pearsonr_values.([hemis{hemi}])(g,:);
    w(idx_not_sig) = 0;
    
    % Make new struct, just containing significant correlation values, and
    % not significant values are zeroed
    sig_mean_pearsonr_values.([hemis{hemi}])(g,:) = w;
    
    end
end

%% Project pearson r values to the surface pre and post fdr correction
mkdir ([output 'images']);
for g = 1:length(age)
    s = surf_subjects{g};
    cd ([data_input s '/surfaces/']);
    surf_rs = SurfStatReadSurf({'lh.smoothwm.native.rsl.obj','rh.smoothwm.native.rsl.obj'});
   %Pre-FDR correction
    pearsonr = [mean_pearsonr_values.l(g,:) mean_pearsonr_values.r(g,:)];
    figure; SurfStatView(pearsonr,surf_rs);
    colormap(RWB)
    SurfStatColLim([-1, 1])
    exportgraphics(gcf, strcat([output 'images/' num2str(age(g))...
        '_' layer{l} '_' dwi_metric{dwi} '_vs_' surf_metric{sm} '_' num2str(24+week) 'w_average_corrcoef_prefdr_map.jpg']),'Resolution',300)
    hold off
   %Post-FDR correction
    sigpearsonr = [sig_mean_pearsonr_values.l(g,:) sig_mean_pearsonr_values.r(g,:)];
    figure; SurfStatView(sigpearsonr,surf_rs);
    colormap(RWB)
    SurfStatColLim([-1, 1])
    exportgraphics(gcf, strcat([output 'images/' num2str(age(g))...
        '_' layer{l} '_' dwi_metric{dwi} '_vs_' surf_metric{sm} '_' num2str(24+week) 'w_average_corrcoef_postfdr_map.jpg']),'Resolution',300)
    hold off
end

