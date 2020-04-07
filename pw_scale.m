function [Q, R] = pw_scale( D, options )
% Scaling method for pairwise comparisons, also for non-balanced
% (incomplete) designs.
%
% [Q, R] = pw_scale( D, options )
%
% D - NxN matrix with positive integers. D(i,j) = k means that the
%     condition i was better than j in k number of trials.
% options - a cell array with the options. Currently recognized options:
%	   'prior' - type of the distance prior in the available options are:
%
%                'none': do not use prior;
%                'bounded': unsurance that the distance between quality 
%                scores is within a bounded range, adaptively selected in
%                each iteration of the MLE optimization;
%                'gaussian': the normalised sum of probabilities of 
%                observing a difference for all compared pairs of conditions.
%
%                Set to 'bounded' by default. Bounded is faster to compute
%                than Gaussian and was marginally worse in the simulation
%                results.
%
%       'regularization' - Since the quality scores in pairwise comparisons
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
% Q - the JOD value for each method. The difference of 1 corresponds to
%     75% of answers selecting one condition over another.
% R - matrix of residuals. The residuals are due to projecting the
%     differences into 1D space. These are the differences between the
%     distances corresponding to the probabilities observed in the
%     experiment and the resulting distances after scaling. 
%
% The condition with index 1 (the first row in D) has the score fixed at value
% 0. Always put "reference" condition at index 1. The measurement error is
% the smallest for the conditions that are the closest to the reference
% condition.
%
% The method scaled the data by solving for maximum-likelihood-estimator
% explaining the collected data. The principle is similar to Bradley-Terry
% model, but Gaussian and not logistic function is used to model relation
% between probabilities and distances in the scaled space. Similar
% scaling was proposed in:
%
% Silverstein, D., & Farrell, J. (2001). Efficient method for paired
% comparison. Journal of Electronic Imaging, 10, 394-398. doi:10.1117/1.1344187
%
% However, the method contains a number of extensions improving robustness,
% eliminating bias and allowing computing confidence intervals (with
% pw_scale_bootstrp function).
%
% Author: Rafal Mantiuk

% Revision history
% 2016-03-19 - Fixed graph connectivity patch; Replaced UA weights with
%              a conditional prior
% 2017-09-13 - Refined the prior and code simplification

% All elements must be non-negative integers 
if any(isinf(D(:))) || any(floor(D(:)) ~= D(:)) || any(D(:)<0)
    error( 'Matrix of comparisons contains invalid inputs');
end

if( ~exist( 'options', 'var' ) )
    options = {};
end

opt = struct();

% We don't the prior by default
opt.prior = 'bounded';
opt.regularization='mean0';
for kk=1:2:length(options)
    if( ~isfield( opt, options{kk} ) )
        error( 'Unknown option %s', options{kk} );
    end
    
    switch options{kk}
        case 'prior'
            if ~ismember( options{kk+1}, { 'none', 'bounded', 'gassian' } )
                error( 'The "prior" option must be "none", "bounded", or "gassian"' );
            end
        case 'regularization'
            if ~ismember( options{kk+1}, { 'mean0', 'fix0' } )
                error( 'The "regularization" option must be "mean0" or "fix0"' );
            end
    end
    
    opt.(options{kk}) = options{kk+1};
    
end
 

if( size(D,1) ~= size(D,2) )
    error( 'The comparison matrix must be square' );
end

% The number of compared conditions
N = size( D, 1 ); 

% Change votes into probabilities - also for incomplete design
M = D./( D + D');
M(isnan(M)) = 0.5;

% inverse cummative distrib, from ISO 20462
Q = sum( -(12/pi * asin( sqrt(M) ) - 3)/2, 1 )';

% find unanimous (UA) and non-unanimous relations (NUA) and build a graph
NUA = (D>0) .* (D'>0); 
G = (D+D')>0;
UA = G - NUA;

% find connected components
group = 0;
node_gr = zeros(1,N);
for rr=1:N    
    if( node_gr(rr) == 0 )
        group = group + 1;
    end
    node_gr = connected_comp( G, node_gr, rr, group );
end

% add links between disconnected components
%[tmp, ord] = sort( Q, 'ascend' );
Ng = max(node_gr); % how many disconnected components
if( Ng > 1 )
    warning( 'There are %d disconnected components in the comparison graph. Some quality scores cannot be accurately computed', Ng );
    % Find the highest quality condition in each disconnected component
    Cb = zeros( Ng, 1 );
    for kk=1:Ng
        Cb(kk) = find( Q == max(Q(node_gr==kk)) & (node_gr==kk)', 1 );
    end

    % Link all highest quality conditions and make them equivalent
    for kk=1:Ng
        for jj=(kk+1):Ng
            D(Cb(kk),Cb(jj)) = 1;
            D(Cb(jj),Cb(kk)) = 1;
            NUA(Cb(kk),Cb(jj)) = 1;
            NUA(Cb(jj),Cb(kk)) = 1;
        end
    end
end

options = optimset( 'Display', 'off', 'LargeScale', 'off' );

D_sum = D + D';
Dt = D';
nnz_d = (D_sum)>0;

% number of pairs compared at least once
comp_made = sum(nnz_d(:));

% Comparison matrix where we shift unanimous answers to the closest
% non-unanimous solution
D_wUA = D;
% Shift anonimous answers equal to 0 to 1
D_wUA(UA==1 & D==0) = 1;
% Substract 1 from the rest of anonimous answers
D_wUA(UA==1 & D~=0) = D_wUA(UA==1 & D~=0) - 1;

f = @(x)exp_prob(x);


if strcmp( opt.regularization, 'mean0' )
    % The methods tend to be more robust if starting scores are 0
    Q_0 = zeros(N,1);
else
    Q_0 = zeros(N-1,1);
end


Q = fminunc( f, Q_0, options );

if ~strcmp( opt.regularization, 'mean0' )
    % Add missing leading 0-score for the first condition (not optimized)
    Q = cat( 1, 0, Q );
end

% Calculate the matrix of residuals
JOD_dist_fit = repmat( Q, [1 N] ) - repmat( Q', [N 1] ); % Compute the distances
JOD_dist_data = norminv( D./D_sum, 0, 1.4826 );

R = NaN( size(D) );
valid = nnz_d & NUA;
R(valid) = JOD_dist_fit(valid) -  JOD_dist_data(valid);

    function P = exp_prob( q_trunc )
        if strcmp( opt.regularization, 'mean0' )
            q = q_trunc;
        else
            q = cat( 1, 0, q_trunc ); % Add the condition with index 1, which is fixed to 0
        end
                        
        sigma_cdf = 1.4826; % for this sigma normal cummulative distrib is 0.75 @ 1
        Dd = repmat( q, [1 N] ) - repmat( q', [N 1] ); % Compute the distances
        Pd = normcdf( Dd, 0, sigma_cdf ); % and probabilities  

        % Compute likelihoods
        prob = Pd(nnz_d);
        p = prob.^D(nnz_d) .*(1-prob).^Dt(nnz_d);        

        L_reg = 0; % regularization loss term
        if strcmp( opt.regularization, 'mean0' )
            L_reg = 0.01 * sum(mean(q).^2);
        end 
        
        % Compute prior
        switch opt.prior
            case 'gaussian'
                
                prior = zeros(comp_made,1);
                for zz=1:N
                    for hh=1:N

                        n = D_sum(zz,hh);
                        %If the comparison has been performed
                        if n>0
                            k = D_wUA(zz,hh);
                            % Compute the probability of each distance
                            % according to all our answers
                            aux = prob.^k .* (1-prob).^(n-k);
                            prior = prior + (aux/sum(aux));

                        end
                    end
                end
                % The mean likelihood per answer is our prior (i.e., we compute
                % the probability of observing a certain distance according to
                % the rest of the answers in our comparison matrix)
            
            case 'bounded'
                q_range = max(q)-min(q);
                n_e = q_range+1;
                prior = max( NUA(nnz_d), 1/n_e - abs(D(nnz_d))/n_e.^2 );
                
            case 'none'
                prior = ones(comp_made,1);
            otherwise
                error( 'Unknown prior option %s', opt.prior );
        end
        
        P = -sum( log( max( p, 1e-400)  ) + ...
            log( max( prior + 0.1, 1e-400) ) ) + ...
            L_reg;

    end

function node_gr = connected_comp( G, node_gr, node, group )
if( node_gr(node) ~= 0 )
    return;
end
node_gr(node) = group;
for nn=find(G(node,:))
    node_gr = connected_comp( G, node_gr, nn, group );
end
end

end
