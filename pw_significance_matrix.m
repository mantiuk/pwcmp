function [p, D] = pw_significance_matrix( jod, pw_stats )
% Get the p-values for the two-tailed z-test between all pairs of
% conditions.
%
% See Section 7.2 "Statistical difference between two conditions" in
% Perez-Ortiz, Maria, and Rafal K. Mantiuk. 
% “A Practical Guide and Software for Analysing Pairwise Comparison Experiments.” 
% ArXiv Preprint, December 11, 2017. http://arxiv.org/abs/1712.03686.
%
% Example:
% [R, Rs] = pw_scale_table( D, [], { 'condition_1', 'condition_2' }, 'user_id', 'selection' );
% [p, D] = pw_significance_matrix( R.jod, Rs{1}.stats ); 

c = pw_stats.jod_cov;
N = size(c,1);

% Variance of the differences between all pairs
V = repmat( diag(c), [1 N] ) + repmat( diag(c)', [N 1] ) - 2*c;

% Diffences between all pairs
D = repmat( jod, [1 N] ) - repmat( jod', [N 1] );

p = normcdf( -abs(D), 0, sqrt(V) );

end