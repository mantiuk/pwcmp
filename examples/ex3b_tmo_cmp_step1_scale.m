% The same as ex3a_tmo_cmp_step1_scale.m, but using utility function
% pw_scale_table()

if( ~exist( 'pw_scale', 'file' ) )
    addpath( fullfile( pwd, '..' ) );
end

D = readtable( 'ex3_tmo_cmp_data.csv' );
D = D( ~strcmp( D.scene, 'peniches' ), : );

% If there are mutiple columns describing a condition (or a group), you can
% concatenate them into a single column, for example:
%
% D.condition_A = strcat( D.distortion_A, '_', D.level_A );
% D.condition_B = strcat( D.distortion_B, '_', D.level_B );

[R, Rs] = pw_scale_table(D, 'scene', { 'condition_1', 'condition_2' }, ...
    'observer', 'selection', 'bootstrap_samples', 500, 'do_all', true );

R = renamevars(R, "condition", "tmo");

C = unique(R.tmo);
writetable(R, 'ex3_tmo_cmp_scaled.csv');
save( 'ex3_tmo_cmp_scaled.mat', 'Rs', 'C' );
