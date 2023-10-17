# DWI2SURF Analysis Pipeline: 
## Software: Matlab2023a
 
# Requirements.
1.     Input data: dHCP-fetal-dwi-and-surf/sub-CC****_ses-****/
/surfaces/: White matter surface reconstruction files, registered to central template (vertex correspondence between all subjects) then warped back to native T2 space, and converted into cartesian coordinates. See Surface Processing Steps for info on how to generate them. (‘${file}/surfaces/?h.smoothwm.native.rsl.coord.obj’)
/volumes/: Co-registration between DWI and T2 (DWI volumes must be in native T2 space)(‘${file}/volumes/tissue_to_t2.nii’)
2.     Sulcal labels of major primary sulci from 28w template.  
3.     List of subject surfaces to use as example for each gestational week (to make figures) (‘example_surface_subjects.csv’).
4.     Analysis output sub-folders within ‘dHCP-fetal-final-cohort’
 
# Step-wise scripts:
## Extract diffusion values beneath the surface, in the cortical plate and subplate. 
 
Master_dwi_to_text.m – extracts values from diffusion maps (tissue/fluid/fa/md), projecting inside and outside surface boundary. (functions- value_under_surface.m, value_outer_surface.m)
 
Output: all_dwi_metric_x2.mat
 
## Extract surface metrics across all vertices and save in struct to match formatting of diffusion values.
 
Master_surf_metrics.m – extracts values from surface reconstructions (depth/curvature/area)
 
Output: all_surf_metric.mat
 
## Calculate the weekly average values in the diffusion and surface metrics
 
Master_weekly_means.m -calculates mean of surface and dwi metrics in each gestational week & projects to surface
(functions- colourscale.m, makes colour bars to use for plotting beforehand)
 
Output: weekly_mean_maps
 
## Age-related trends in diffusion and surface metrics, pearson r correlation values on surface.
 
Master_dwi_and_surf_vs_age.m -pearson r value between age and dwi/surf metrics plotted and projected to the surface
 
Output: dwi_vs_age_corr_maps, surf_vs_age_corr_maps

## Neighbourhood analysis – local relationships (coupling) between diffusion and surface metrics in individual subjects.
 
Master_neighbourhood.m – finds neighbours for the surface vertex mesh
Master_coupling.m – uses the neighbours’ output from 5.1 and correlates values within subject.
 
### Functions:
FDR.m, RWB.mat
Output: 
Coupling/images
 
## Age mismatch coupling – neighbourhood analysis for older sulcal depth vs. younger tissue fraction
 
Master_age_mistmatch.m – takes mean surface metrics in older gestational week and correlates them with younger native diffusion values. 
 
### Functions:
Colourscale.m
Output: 
age-mismatch

## Whole-brain age mismatch coupling- a modified script derived from master_coupling.m that uses average sulcal depth in older gestational weeks, correlates it against native DWI metrics
 
Whole_brain_coupling_mismatch.m
Output:
Coupling-mismatch-whole-brain
