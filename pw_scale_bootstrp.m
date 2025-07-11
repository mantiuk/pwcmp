function [jod, stats] = pw_scale_bootstrp( MM, bootstrap_samples, opt )
arguments
    MM {mustBeMatrix}
    bootstrap_samples {mustBeNonnegative, mustBeInteger} = 1
    opt.display {mustBeMember(opt.display, {'info', 'none'})} = 'info'
    opt.alpha {mustBeInRange(opt.alpha, 0, 1)} = 0.05
    opt.use_parallel {mustBeMember(opt.use_parallel, {'always', 'never'})} = 'always'
    opt.regularization {mustBeMember(opt.regularization, {'mean0', 'fix0'})} = 'mean0'
    opt.prior {mustBeMember(opt.prior, {'gaussian', 'none'})} = 'gaussian' 
end
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
%      'display' - set to 'none' for quiet operation, and 'info' to show
%                some extra information. 'info' is the default option.
%      'alpha' - the 'alpha' value for condidence interval. Default value
%                of 0.05 results in 95% confidence intervals.
%      'use_parallel' - use parallel processing for bootstrapping. Possible
%                values: 'always' (default) or 'never'.
%
%	   'prior' - type of the distance prior in the available options are:
%                'none': do not use prior;
%                'gaussian': the normalised sum of probabilities of 
%                observing a difference for all compared pairs of conditions.
%
%                Set to 'gaussian' by default. 
%      
%      'regularization' - Since the quality scores in pairwise comparisons
%                are relative and the absolute value cannot be obtained, it
%                is necessary to make an assumption how the absolute values
%                are fixed in the optimization. The two options are:
%
%                'mean0' - add a regularization term that makes the mean
%                JOD value equal to 0. 
%                'fix0' - fix the score of the first condition to 0. That
%                score is not optimized. 
%
%                The default is 'mean0'. 'mean0' results in a reduced
%                overall estimation error as compared to 'fix0'. 'fix0' is 
%                useful when one of the conditions is considered a
%                baseline or a reference. The conditioons closer to that
%                reference will be estimated with higher accuracy. 
%
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
% stats.bstrp - the complete set of bootstrap samples as a [B,C] matrix,
%      where B is the number of bootstrat samples and C is the number of
%      conditions


if( bootstrap_samples==0 )
    bootstrap_samples = 1;
end

opt.display_level = strcmp( opt.display, 'info' );
    
N = sqrt( size(MM,2) );

M = reshape( sum(MM,1), N, N );

M = round_simulated_ties(M);

[jod, R] = pw_scale( M, {'prior', opt.prior, 'regularization', opt.regularization});
Rv = abs(triu(R,1));
if( opt.display_level > 0 )
    fprintf( 1, 'Residual due to scaling: mean = %g; min = %g; max = %g\n', mean(Rv(:), "omitmissing"), min(Rv(:)), max(Rv(:)) );
end

stats = struct();

if( bootstrap_samples == 1 || size(MM,1) < 2 )
    stats.jod_low = jod;
    stats.jod_high = jod;
    stats.jod_cov = zeros(N,N);
    return;
end

if( opt.display_level > 0 )
    fprintf( 1, 'Generating %d bootstrap samples. It can take a while...\n', bootstrap_samples );
end

% Use if parallel proc toolbox available
options = statset( 'UseParallel', opt.use_parallel ); 
%options = statset(  );
bstat = bootstrp( bootstrap_samples, @boot_jod, MM, 'Options', options );

stats.bstrp = bstat;

% Test if each JOD-scaled point belongs to standard distribution
H_p = 0;
if strcmp( opt.regularization, 'fix0' )
    first_cond = 2; % Skip the anchor condition, as std=0
else
    first_cond = 1;
end

for kk=first_cond:size(bstat,2)
    xb = bstat(:,kk);
    x = (xb-mean(xb))/std(xb);
    H_p = H_p + (1-kstest( x ));
end
if( opt.display_level > 0 )
    fprintf( 1, '%d out of %d JOD-points have a standard normal distribution (Kolmogorov-Smirnov test)\n', H_p, size(bstat,2)-first_cond+1 );
end

stats.jod_low = prctile( bstat, opt.alpha*100/2 )';
stats.jod_high = prctile( bstat, 100 - opt.alpha*100/2 )';
stats.jod_cov = cov( bstat )';

    function Q = boot_jod( MM_bst )
        M = reshape( sum(MM_bst,1), N, N );
        M = round_simulated_ties( M );
        Q = pw_scale( M, {'prior', opt.prior, 'regularization', opt.regularization} );
    end

end

function M = round_simulated_ties(M)
    % In simmulation we allow ties i.e. an option of not giving preference
    % to one of two conditions results in 0.5 is assigned to both in a given
    % round. Since fractions are not allowed in the scaling, rounding is 
    % performed. The function rounds non-integer entries of the pwc matrix 
    % to the smaller or larger integer (randomly) but agreeing with the 
    % number of comparisons.
    N = size(M,1);
    aux = triu((randi(2,[N,N])-1),1);
    M = (aux + triu(-(aux-1),1)').*round(M) + ~(aux + triu(-(aux-1),1)').*floor(M);
end
