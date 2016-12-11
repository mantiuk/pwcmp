function [p, D] = pw_significance_matrix( jnd, pw_stats )

c = pw_stats.jnd_cov;
N = size(c,1);

% Variance of the differences between all pairs
V = repmat( diag(c), [1 N] ) + repmat( diag(c)', [N 1] ) - 2*c;

% Diffences between all pairs
D = repmat( jnd, [1 N] ) - repmat( jnd', [N 1] );

p = normcdf( -abs(D), 0, sqrt(V) );

end