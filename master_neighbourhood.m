%% FINAL VERSION of script to find neighbours
% Requirements:
% mris_convert -v lh.smoothwm.native.rsl.asc
% lh.smoothwm.native.rsl.neighbours.asc (for left and right)

% mris_convert writes out neighbors of a vertex in each row. 
% The 1st column is the vertex number, the 2nd col is the number of neighbors,
% the remaining cols are the vertex numbers of the neighbors.  
% Note: there can be a different number of neighbors for each vertex.

% Convert the neighbours.asc to excel file and fill empty cells so matlab can read it as a table

% This script will differentiate between neighbors of various distances
% so that weights can be given to neighbors of different degrees
% each set of neighbors is neighbors of only that degree

% Set the directory containing the 'template' data i.e
% .neigbours.csv files
main = ('/Users/sianwilson/Desktop/BCH/MASTER/');

% Set paths to directories
templatedir = ([main 'dHCP-fetal-final-cohort/coupling/']);
mkdir(templatedir);
surf_dir = ([main 'dHCP-fetal-final-cohort/dHCP-fetal-dwi-and-surf/']);
addpath([main 'Scripts']);

% Define the hemisphere names
hemis = {'l', 'r'};

% Define max degree of separation (check what was used in neighbourhood.m
% script)
k = 5;

% Initialise struct with cells to keep info
N_info.lh = cell(length(k));
N_info.rh = cell(length(k));

%Define the subject
subject = ('sub-CC01226XX14_ses-150730');
   
%Define the path to subject surface data
surf_path = ([surf_dir,subject,'/surfaces/']);

%% Find first neighbours
% Loop through each hemisphere
for h = 1:length(hemis)
    hemi = hemis{h};
    N_info.(hemi) = cell(length(k));
    file = SurfStatReadSurf([surf_path char(hemi) 'h.smoothwm.to31.rsl.obj']);
    N_info.(hemi){1} = cell(length(file.coord),1);
    for V = 1:length(file.coord)
        N_info.(hemi){1}{V} = find_neighbour_info(V,file);
    end
end
save(fullfile(templatedir, ['lh_rh.', num2str(k), 'neighbors_1st_degree.mat']), 'N_info');    

%% Find neighbours of neighbours up to k degrees of separation
vertices = num2cell(1:length(file.coord))';

for h = 1:length(hemis)
    hemi = hemis{h};
    fprintf (hemi)
    for deg = 2:k
        fprintf('%d\t', deg);
        N_info.(hemi){deg} = cell(length(file.coord),1);
        temp_array = cat(2,N_info.(hemi){:},vertices);
       for V = 1:length(file.coord)
            temp = arrayfun(@(x) find_neighbour_info(x,file), N_info.(hemi){deg-1}{V},'UniformOutput', false);
            A = unique(cell2mat(temp));
            B = cat(1,temp_array{V,:});
            uni = setdiff(A, B);
            N_info.(hemi){deg}{V} = [N_info.(hemi){deg}{V} ;uni]; 
        end
    end
    save(fullfile(templatedir, ['lh_rh.', num2str(k), 'neighbors_of_neighbours.mat']), 'N_info');    
end

%% Reformat the neighbours.mat file
neighbours = cell(1,k); 

for hemi = 1:length(hemis)
    for i = 1:k
    neighbours{i} = cellfun(@double, N_info.([hemis{hemi}]){i},'UniformOutput', false);
    end
    neighbours_cat.([hemis{hemi}]) = arrayfun(@(row) horzcat(neighbours{row,:}), 1:size(neighbours,1), 'UniformOutput', false);
    neighbours_cat.([hemis{hemi}]) = neighbours_cat.([hemis{hemi}]){1};
    save(fullfile(templatedir, ['cat_format_lh_rh.', num2str(k), 'neighbors_of_neighbours.mat']), 'neighbours_cat');
end

    