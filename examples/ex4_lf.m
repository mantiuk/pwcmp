% This is an example of scaling a large-size experiment (26580 comparisons) 
%
% Set boostrap_samples = 0 to compute JODs but no confidence intervals. 
%
% The data for this example comes from the paper:
% Adhikarla, V. K., Vinkler, M., Sumin, D., Mantiuk, R. K., Myszkowski, K., Seidel, H.-P., & Didyk, P. (2017). 
% Towards a quality metric for dense light fields. 
% Computer Vision and Pattern Recognition (CVPR), 58–67.
%
% Project web page: http://lightfields.mpi-inf.mpg.de/

if( ~exist( 'pw_scale', 'file' ) )
    addpath( '../' );
end

D = readtable( 'lf_quality_pwcmps.csv' );

% A condition is a combination of distotion type and distortion level. 
% Create a column for each compared condition.
D.condition_A = strcat( D.dist_type1, '_', num2str( D.dist_level1 ) );
D.condition_B = strcat( D.dist_type2, '_', num2str( D.dist_level2 ) );
D.is_A_selected = (D.selected==1);

[R, Rs] = pw_scale_table(D, 'scene', { 'condition_A', 'condition_B' }, ...
    'observer', 'is_A_selected', 'bootstrap_samples', 500, 'do_all', true, 'regularization', 'mean0', 'anchor_condition', 'Reference_ 0', 'use_parallel', 'never' );

writetable( R, 'lf_quality_scaled.csv' );

