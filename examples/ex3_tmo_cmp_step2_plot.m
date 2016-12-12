% Part 2 of the ex5. Run ex5_tmo_cmp_scale.m first to generate scale
% results
%
% This example shows two ways to plot data

% Method 1. As a bar plot with error bars

figure(1)

D = dataset( 'file', 'ex3_tmo_cmp_scaled.csv', 'delimiter', ',' );

D.tmo = strrep( D.tmo, '_', '-' );
SCs = unique( D.scene );
C = unique( D.tmo );


% Offset all results so that they are all positive starting from 1
jod_offset = 1 - min(D.jod);

N = length(C); % number of TMOs (conditions)
M = length(SCs); % number of scenes

JODs = zeros(M, N);
error_bars = zeros(M, N, 2);

for kk=1:length(SCs)
    for ll=1:length(C)
        ind = find( strcmp(D.scene,SCs{kk}) & strcmp(D.tmo,C{ll}) );
        JODs(kk,ll) = D.jod( ind ) + jod_offset;
        error_bars(kk,ll,1) = D.jod_high( ind )-D.jod( ind );
        error_bars(kk,ll,2) = D.jod( ind )-D.jod_low( ind );
    end
end
barwitherr( error_bars, JODs );
set( gca, 'XTickLabel', SCs );
set( gca, 'YTick', 1:10 );
legend( C, 'Location', 'best' );
grid on


% Method 2. As "triangle" plots showing whether the differences between
% conditions are statistically significnat. Red-dashed line - not
% significant; blue-continous line - significant.

figure(2);

res = load( 'ex3_tmo_cmp_scaled.mat' );
M = length(res.Rs); % Number of scenes

clf
for kk=1:M
    subplot( ceil(M/2), 2, kk );
    pw_plot_ranking_triangles( res.Rs{kk}.jod, res.Rs{kk}.stats, res.C );    
    title( res.Rs{kk}.scene );
    xlabel( 'JODs' );
end
