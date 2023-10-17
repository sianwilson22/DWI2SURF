%% This script converts the diffusion metric maps (in native t2 volume space) to a txt file, 
% by interpolating from the surface normal to a ribbon of voxels beneath
% the white matter surface
% The output is two txt files, representing values in left and right hemispheres

% User must input the list of subject ids, metric of interest (tissue/fluid/fa/md), and the
% corresponding file name (fa_to_t2.nii or fluid_to_t2.nii)

% Input files are structured so that 'dir' contains subject folders eg.
% subject/volumes/tissue_to_t2.nii 
% subject/surfaces/lh.smoothwm.native.rsl.coord.obj

%Output is images of diffusion metrics projected to the surface, saved in
%img_outdir and .txt files saved in txt_outdir. 
%For future analysis, a struct 'all_dwi_metric.mat' is constructed- the
%order is determined by variable: dwi_metrics = {'tissue';'fluid';'fa';'md'};

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

% Set path for output directory
date = today("datetime");
str_date = datestr(date);

mkdir([main_input 'dwi-inner-outer-to-surf-txt-' str_date]);
txt_outdir = ([main_input 'dwi-inner-outer-to-surf-txt-' str_date '/']);

%% Diffusion values beneath the white matter surface
%Define metric being used
dwi_metric = {'tissue';'fluid';'fa';'md'};
hemis = {'l','r'};
data_input = [main_input 'dHCP-fetal-dwi-and-surf/'];

for hemi = 1:2
    all_dwi_metric.inner.one.(hemis{hemi}) = cell(length(dwi_metric),1);
    all_dwi_metric.inner.two.(hemis{hemi}) = cell(length(dwi_metric),1);
    all_dwi_metric.outer.one.(hemis{hemi}) = cell(length(dwi_metric),1);
    all_dwi_metric.outer.two.(hemis{hemi}) = cell(length(dwi_metric),1);

    for dwi = 1:4
    all_dwi_metric.inner.one.(hemis{hemi}){dwi} = zeros(length(no_space_list),vertices);
    all_dwi_metric.inner.two.(hemis{hemi}){dwi} = zeros(length(no_space_list),vertices);
    all_dwi_metric.inner.one.(hemis{hemi}){dwi} = zeros(length(no_space_list),vertices);
    all_dwi_metric.outer.two.(hemis{hemi}){dwi} = zeros(length(no_space_list),vertices);

        for n = 1:length(no_space_list)
            %Define the subject
            x = string(no_space_list{n});
            subject = regexprep(x, '\s', '' );
            
            %Define the path to subject surface data
            surf_path = strcat(data_input,subject,"/surfaces/");
            cd (surf_path)
            
            %Read in surface data
            s =  SurfStatReadSurf1([hemis{hemi} 'h.smoothwm.native.rsl.coord.obj']);
            
            %Define path to volume data
            vol_path = strcat(data_input,subject,"/volumes/");
            cd (vol_path)
            dwi_vol = [dwi_metric{dwi} '_to_t2.nii'];
            %gunzip ([dwi_vol '.gz'])
            vol = SurfStatReadVol1(dwi_vol); %CHANGE HERE
            map = vol.data;
            
            %Calculate diffusion values on left & right surface using function:
            %value_under_surface.m
            [val_inner1,val_inner2] = value_under_surface(s, map);
            [val_outer1,val_outer2] = value_outer_surface(s, map);
            
            %Fill in the struct for future use
            all_dwi_metric.inner.one.(hemis{hemi}){dwi}(n,:) = val_inner1;
            all_dwi_metric.inner.two.(hemis{hemi}){dwi}(n,:) = val_inner2;

            all_dwi_metric.outer.one.(hemis{hemi}){dwi}(n,:) = val_outer1;
            all_dwi_metric.outer.two.(hemis{hemi}){dwi}(n,:) = val_outer2;

            %Save the inner white matter and outer values for each subject
            cd(txt_outdir)
            filename.inner1 = strcat(subject,"_",dwi_metric{dwi},".", hemis{hemi},"h.smoothwm.native.inner.1.diffusion");
            filename.inner2 = strcat(subject,"_",dwi_metric{dwi},".", hemis{hemi},"h.smoothwm.native.inner.2.diffusion");
            filename.outer1 = strcat(subject,"_",dwi_metric{dwi},".", hemis{hemi},"h.smoothwm.native.outer.1.diffusion");
            filename.outer2 = strcat(subject,"_",dwi_metric{dwi},".", hemis{hemi},"h.smoothwm.native.outer.2.diffusion"); 
            a = fopen(filename.inner1,'w');
            b = fopen(filename.inner2,'w');
            c = fopen(filename.outer1,'w');
            d = fopen(filename.outer2,'w');
            fprintf(a, '%g \n', val_inner1);
            fprintf(b, '%g \n', val_inner2);
            fprintf(c, '%g \n', val_outer1);
            fprintf(d, '%g \n', val_outer2);
            fclose(a);
            fclose(b);
            fclose(c);
            fclose(d);
        end
    end
end
cd (data_input);    
save('all_dwi_metric_x2.mat', 'all_dwi_metric');