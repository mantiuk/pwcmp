function L = pw_outlier_analysis( M_list )
% A simple outlier analysis for a pairwise comparison data
%
% L = pw_outlier_analysis( MM )
% L = pw_outlier_analysis( { MM1, MM2, ..., MMn } )
%
% The function takes as input matrix MM, which contains responses for each
% observer in seperate rows. The format of that matrix is the same as for
% pw_scale_boostrp function.
%
% The function can also take as input a list of MM matrices. This is useful
% if a pairwise comparison experiment was performed for several scenes but
% with the same set of observers. If the list is provided, likelihoods are
% summed across all scenes, give more reliable data. 
% 
% The function returns a vector of log-10 likelihoods, one entry for each
% observer. The values represent a log-10 likelihood of observing the data 
% for k-th observer (participant) if the true JOD values are given by scaling
% data for all other participants. 
% 
% The function also displays log-10 likelihood and log-ratio difference
% for each observer. Large positive log-ratio difference (e.g. greater than
% 1) could indicate that a particular observer is an outlier. 
% 
% The result of this analysis should be interpreted with care: small
% differences in log-10 likelihood normally occur due to random effects.
% Observers can be rejected only if the (negatiuve) log-likelihood is much
% smaller then for other observers.
%
% Author: Rafal Mantiuk

if( ~iscell( M_list ) )
    M_list = { M_list };
end

N = sqrt(size(M_list{1},2));  % The number of conditions
assert( (N - floor(N)) == 0 );
N_obs = size(M_list{1},1);  % The number of observers

L = zeros(N_obs,1); % Log-likelihood per observer

for ss=1:length(M_list)
    
    MM = M_list{ss};
    assert( size(MM,1) == N_obs ); % The number of observers must be the same for each scene
    
    % For each observer
    for oo=1:N_obs
        
        sel = true(size(MM,1),1);
        sel(oo) = false;
        MMp = MM(sel,:);
        MMref = reshape( sum(MMp), [N N] );
        
        %    [jod, stats] = pw_scale_bootstrp( MMp, bootstrap_samples, options );
        jod = pw_scale( MMref );
        
        % Given the posterior for all other observers, how likely is that the
        % observer "oo" result was generated from the same distribution
        
        M = reshape( MM(oo,:), [N N] );
        
        sigma_cdf = 1.4826; % for this sigma normal cummulative distrib is 0.75 @ 1
        Dd = repmat( jod, [1 N] ) - repmat( jod', [N 1] ); % Compute the distances
        Pd = normcdf( Dd, 0, sigma_cdf ); % and probabilities
        
        %Dt = D';
        
        D_sum = M + M';
        
        % Compute likelihoods for N<=30 and N>30
        p = binopdf( M, D_sum, Pd );
        
        % Which probabilities to consider
        sel_p = (triu(ones(size(p)),1)==1) & (D_sum>0);
        p_all = p(sel_p);
        
        L(oo) = L(oo) + sum( log10( max( p_all, 1e-200) ) );
                
    end
    
end

for oo=1:N_obs
    sel = true(N_obs,1);
    sel(oo) = false;
    L_median = median(L(sel));
    lld = log(L(oo)/L_median);
    display( sprintf( 'Log-10 likelihood for observer %d = %g;  log-ratio difference to mean = %g', oo, L(oo), lld ) );
end

end