%Directory containing the subject surface files
main = ('/Users/sianwilson/Desktop/BCH/MASTER/');

%Add path to functions
addpath([main 'Scripts']);

main_input = ([main 'dHCP-fetal-final-cohort/']);

% Set number of vertices
vertices = 40962;

% Group of subjects
subject_list = importdata ([main 'Subject-Lists/age-sorted-subjects.txt']);
no_space_list = regexprep(subject_list, '\s', '');

%Create struct containing all surface metrics for each subject
%Define metric being used
surf_metric = {'depth';'area';'mc'};
hemis = {'l','r'};
data_input = [main_input 'dHCP-fetal-dwi-and-surf/'];

for hemi = 1:2
    all_surf_metric.(hemis{hemi}) = cell(length(surf_metric),1);
 
    for surf = 1:3
    all_surf_metric.(hemis{hemi}){surf} = zeros(length(no_space_list),vertices);

        for n = 1:length(no_space_list)
            %Define the subject
            x = string(no_space_list{n});
            subject = regexprep(x, '\s', '' );
            
            %Define the path to subject surface data
            surf_path = strcat(data_input,subject,"/surfaces/");
            cd (surf_path)
            
            %Load file
            a = load([ hemis{hemi} 'h.' surf_metric{surf}]);
            %Fill in the struct for future use
            all_surf_metric.(hemis{hemi}){surf}(n,:) = a;
        end
    end
end
cd (data_input);    
save('all_surf_metric.mat', 'all_surf_metric');