function [Q, R] = pw_scale( D )
% Scaling method for pairwise comparisons, also for non-balanced
% (incomplete) designs.
%
% [Q, R] = pw_scale( D )
%
% D - NxN matrix with positive integers. D(i,j) = k means that the
%     condition i was better than j in k number of trials.
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


if( size(D,1) ~= size(D,2) )
    error( 'The comparison matrix must be square' );
end

N = size( D, 1 );  % The number of compared conditions

% Change votes into probabilities - also for incomplete design
M = D;
for rr=1:N
    for cc=(rr+1):N
        C = M(rr,cc)+M(cc,rr);
        if( C == 0 )
            M(rr,cc) = 0.5;
            M(cc,rr) = 0.5;
        else
            M(rr,cc) = M(rr,cc)/C;
            M(cc,rr) = M(cc,rr)/C;
        end
    end
end
M(eye(N)==1) = 0.5;

% inverse cummative distrib, from ISO 20462
S = 12/pi * asin( sqrt(M) ) - 3;

Q = sum( -S/2, 1 )';

% find non-unanimous relations and build a graph
NUA = zeros(N,N);
G = zeros(N,N); %graph
for rr=1:N
    for cc=1:N
        if( D(rr,cc)>0 && D(cc,rr)>0 )
            NUA(rr,cc) = 1;
        end
        if( D(rr,cc)>0 || D(cc,rr)>0 )
            G(rr,cc) = 1;
        end
    end
end

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
nnz_d = (D_sum)>0;

% Use binomial distribution when N<=30, Gaussian otherwise (is faster and is a good
% approximation)
nnz_bino = (D_sum<=30) & nnz_d;
nnz_gauss = (D_sum>30) & nnz_d;

% Precompute to speed-up computation
NK = zeros(N,N);
for ii=1:N
	for jj=1:N
        if( nnz_bino(ii,jj) )
            NK(ii,jj) = nchoosek( D_sum(jj,ii), D(ii,jj) );
        end
    end
end


% The methods tend to be more robust if starting scores are 0
Q = fminunc( @exp_prob, zeros(N-1,1), options );

% Add missing leading 0-score for the first condition (not optimized)
Q = cat( 1, 0, Q );

% Calculate the matrix of residuals
JOD_dist_fit = repmat( Q, [1 N] ) - repmat( Q', [N 1] ); % Compute the distances
JOD_dist_data = norminv( D./D_sum, 0, 1.4826 );

R = NaN( size(D) );
valid = nnz_d & NUA;
R(valid) = JOD_dist_fit(valid) -  JOD_dist_data(valid);

    function P = exp_prob( q_trunc )
        
        q = cat( 1, 0, q_trunc ); % Add the condition with index 1, which is fixed to 0
                
        q_range = max(q)-min(q);
        
        sigma_cdf = 1.4826; % for this sigma normal cummulative distrib is 0.75 @ 1
        Dd = repmat( q, [1 N] ) - repmat( q', [N 1] ); % Compute the distances
        Pd = normcdf( Dd, 0, sigma_cdf ); % and probabilities
                            
        Dt = D';

        % Compute likelihoods for N<=30 and N>30
        p1 = NK(nnz_bino).*Pd(nnz_bino).^D(nnz_bino).*(1-Pd(nnz_bino)).^Dt(nnz_bino);
        p2 = binopdf( D(nnz_gauss), D_sum(nnz_gauss), Pd(nnz_gauss) );
        
        % The prior 
        n_e = q_range+1;
        p1_prior = max( NUA(nnz_bino), 1/n_e - abs(D(nnz_bino))/n_e.^2 );
        p2_prior = max( NUA(nnz_gauss), 1/n_e - abs(D(nnz_gauss))/n_e.^2 );
        
        P = -sum( log( max( [p1.*p1_prior; p2.*p2_prior], 1e-200) ) );                        
    end


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
