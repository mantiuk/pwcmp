function [jod, stats] = pw_scale_bootstrp( MM, boostrap_samples, options )
% Scaling method for pairwise comparisons with bootstrapped confidence
% intervals.
%
% [jod, stats] = pw_scale_bootstrp( MM )
% [jod, stats] = pw_scale_bootstrp( MM, boostrap_samples )
% [jod, stats] = pw_scale_bootstrp( MM, boostrap_samples, options )
%
% Use this function to scale the data from the pairwise comparison
% experiment in just-objectionable-difference (JOD) units and to find 
% the confidence intervals for JOD values. The function employes 
% bootstrapping and may take longer to complete. It can be accelerated when 
% parralel pool is started (see parpool). 
% 
% MM - KxC matrix with positive integers. Each row of the matrix should
%      contain a comparison matrix for a single observer in 'flattenned' 
%      vectoried format. That is, if a comparison matrix for observer k is
%      N, MM should be initialized: 
%      MM(k,:) = M(:);
%      The number of columns C should be equal to N*N where N is the number
%      of compared conditions. The number of rows K is equal the number of
%      observers. Refer to the documentation of pw_scale on the explanation
%      how matrix M is constricted. 
%
% boostrap_samples - The number of samples to use to compute confidence
%      intervals using bootsrapping. Pass 0, 1 or an empty matrix to
%      disable bootstrapping and compute only scaling. Typical number of
%      bootstrap samples can vary between 500 and 2000 depending on the
%      number of observers. Higher number takes longer to compute but could
%      result in more accurate estimates, especially if the number of
%      observers is high. 
%
% options - a cell array with the options. Currently recognized options:
%      'display' - set to 'none' for quite operation, and 'info' to show
%      some extra information. 'info' is the default option.
%      'alpha' - the 'alpha' value for condidence interval. Default value
%      of 0.05 results in 95% confidence intervals.
%      'use_parallel' - use parallel processing for bootstrapping. Possible
%      values: 'always' (default) or 'never'.
%	   'prior' - Boolean indicating whether to use distance prior in the 
%	   optimization. Set to 'true' by default. 
%
% The function return:
% jod - the JOD assigned to each condition. The firt element will be always
%      0 (refer to pw_scale).
% 
% stats - the structure with additional results:
% stats.jod_high - the upper value of the confidence interval. The height
%      of the bar on the error plot should be stats.jod_high - jod.
% stats.jod_low - the lower value of the confidence interval. The height
%      of the bar on the error plot should be jod - stats.jod_low.
% stats.jod_cov - the covariance matrix for the jod scores. It could be
%      used for statistical testing significant differences.

% Author: Rafal Mantiuk

if( ~exist( 'boostrap_samples', 'var' ) || boostrap_samples==0 )
    boostrap_samples = 1;
end

if( ~exist( 'options', 'var' ) )
    options = {};
end

opt = struct();
opt.display = 'info';
opt.alpha = 0.05;
opt.use_parallel = 'always';
opt.prior = 1;
for kk=1:2:length(options)
    if( ~isfield( opt, options{kk} ) )
        error( 'Unknown option %s', options{kk} );
    end
    opt.(options{kk}) = options{kk+1};
end

if( opt.alpha < 0 || opt.alpha > 1 )
    error( 'The "alpha" parameters must be between 0 and 1' );
end

opt.display_level = strcmp( opt.display, 'info' );
    
N = sqrt( size(MM,2) );

M = reshape( sum(MM,1), N, N );

[jod, R] = pw_scale( round(M), opt.prior );
Rv = abs(R(~isnan(R)));
if( opt.display_level > 0 )
    display( sprintf( 'Residual due to scaling: mean = %g; min = %g; max = %g', mean(Rv), min(Rv), nanmax(Rv) ) );
end

stats = struct();

if( boostrap_samples == 1 || size(MM,1) < 2 )
    stats.jod_low = jod;
    stats.jod_high = jod;
    stats.jod_cov = zeros(N,N);
    return;
end

if( opt.display_level > 0 )
    display( sprintf( 'Generating %d bootstrap samples. It can take a while', boostrap_samples ) );
end

% Use if parallel proc toolbox available
options = statset( 'UseParallel', opt.use_parallel ); 
%options = statset(  );
bstat = bootstrp( boostrap_samples, @boot_jod, MM, 'Options', options );

stats.bstrp = bstat;

% Test if each JOD-scaled point belongs to standard distribution
H_p = 0;
for kk=2:size(bstat,2)
    H_p = H_p + kstest( bstat(kk,:) );
end
if( opt.display_level > 0 )
    display( sprintf( '%d out of %d JOD-points have a standard normal distribution (Kolmogorov-Smirnov test)', H_p, size(bstat,2)-1 ) );
end

stats.jod_low = prctile( bstat, opt.alpha*100/2 )';
stats.jod_high = prctile( bstat, 100 - opt.alpha*100/2 )';
stats.jod_cov = cov( bstat )';

function Q = boot_jod( MM_bst )    
    M = reshape( sum(MM_bst,1), N, N );
    Q = pw_scale( M, opt.prior );
end


end

