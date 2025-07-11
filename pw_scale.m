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
%                'gaussian': the normalised sum of probabilities of 
%                observing a difference for all compared pairs of conditions.
%
%                Set to 'gaussian' by default. 
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
opt.prior = 'gaussian';
opt.regularization = 'mean0';
opt.use_gradients = true; 
for kk=1:2:length(options)
    if( ~isfield( opt, options{kk} ) )
        error( 'Unknown option %s', options{kk} );
    end
    
    switch options{kk}
        case 'prior'
            if ~ismember( options{kk+1}, { 'none', 'bounded', 'gaussian' } )
                error( 'The "prior" option must be "none", "bounded", or "gaussian"' );
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

if strcmp( opt.regularization, 'mean0' )
    % The methods tend to be more robust if starting scores are 0
    Q_0 = zeros(N,1);
else
    Q_0 = zeros(N-1,1);
end

options = optimoptions('fminunc','SpecifyObjectiveGradient', opt.use_gradients, 'Display', 'off');
Q = fminunc( @exp_prob_grad, Q_0, options );

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

    function [P, grad] = exp_prob_grad( q_trunc )
        if strcmp( opt.regularization, 'mean0' )
            q = q_trunc;
        else
            q = cat( 1, 0, q_trunc ); % Add the condition with index 1, which is fixed to 0
        end
                        
        sigma_cdf = 1.4826; % for this sigma normal cummulative distrib is 0.75 @ 1
        Dd = repmat( q, [1 N] ) - repmat( q', [N 1] ); % Compute the distances
        
        % early slicing to save memory
        Dd = Dd(nnz_d);
        [row, col] = find(nnz_d);
        dDd_dq = zeros(comp_made, N);

        dDd_dq(sub2ind(size(dDd_dq), (1:comp_made)', row)) = 1;
        coller = sub2ind(size(dDd_dq), (1:comp_made)', col);
        dDd_dq(coller) = dDd_dq(coller) - 1;

        Pd = normcdf( Dd, 0, sigma_cdf ); % and probabilities

        % Compute likelihoods

        prob = Pd;
        Dn = D(nnz_d);
        Dtn = Dt(nnz_d);
        dprob_dq = normpdf(Dd, 0, sigma_cdf) .* dDd_dq;
        p = prob.^ Dn .*(1-prob).^Dtn;
        %dp_dq = prob.^(max(Dn-1, 0)).*((1-prob).^(max(Dtn-1,0)))...
        %            .*(Dn - Dtn.*prob) .* dprob_dq;
        dp_dq = (Dn.* prob.^(max(Dn-1, 0)) .* ((1-prob).^(max(Dtn,0))) - ...
                    Dtn .* ((1-prob).^(max(Dtn - 1,0))) .* prob.^Dn) .* dprob_dq;

        L_reg = 0; % regularization loss term
        dLreg_dq = 0;
        if strcmp( opt.regularization, 'mean0' )
            L_reg = 0.01 * sum(mean(q).^2);
            dLreg_dq = 0.01 * 2 * mean(q) * q / size(q,1);
        end 
        
        % Compute prior
        switch opt.prior
            case 'gaussian'
                % Vectorize the double for loop                
                pos = find(D_sum > 0);
                n = D_sum(pos)';
                k = D_wUA(pos)';

                epsilon = 1e-8;
                aux = prob.^(k) .* (1-prob).^((n-k));
                sumaux = sum(aux);
                prior = sum((aux+epsilon)./(sumaux+epsilon), 2);

                % part of the derivatives.
                daux_dq_A = (k .* prob.^(max(k-1,0)) .* (1-prob).^((n-k)) + ...
                            prob.^k .* (- n + k) .* (1-prob).^(max(n-k-1,0)));

                part1 = (daux_dq_A*((1+epsilon)./(sumaux+epsilon)')) .* dprob_dq;
                part2 = ((1+epsilon)./(sumaux.^2+epsilon) .*aux) * (daux_dq_A'*dprob_dq);
                dprior_dq = part1 - part2;

                % The mean likelihood per answer is our prior (i.e., we compute
                % the probability of observing a certain distance according to
                % the rest of the answers in our comparison matrix)
            
            case 'bounded'
                error( '"bounded" prior is no longer supported. Use "gaussian" instead.')
%                 q_range = max(q)-min(q);
%                 n_e = q_range+1;
%                 prior = max( NUA(nnz_d), 1/n_e - abs(D(nnz_d))/n_e.^2 );
                
            case 'none'
                prior = ones(comp_made,1);
            otherwise
                error( 'Unknown prior option %s', opt.prior );
        end
        
        dpart_dq = zeros(comp_made, N);
        dpart_dq  = dpart_dq + 1./max( p, 1e-400).*( p > 1e-400)  .* (dp_dq);
        if strcmp( opt.prior, 'gaussian' )
            dpart_dq  = dpart_dq + 1./max(prior + 0.1, 1e-400) .* (prior > -0.1 + 1e-400) .* dprior_dq;
        end

        grad = dLreg_dq - sum(dpart_dq)';
        if ~strcmp( opt.regularization, 'mean0' )
            grad = grad(2:end);
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
