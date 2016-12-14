function [q_s] = pw_simulate_exp( N_exps, q, N_obs, C_mat )
% Monte-Carlo simulation of runnimg 'N_exp' experiment when the true
% quality scores are provided in the vector 'q'.
%
% [q_s] = pw_simulate_exp( N_exps, q, N_obs )
% [q_s] = pw_simulate_exp( N_exps, q, N_obs, C_mat )
%
% N_exp - how many experiments to simulate (>=1000 should work)
% q - a row vector with the true value of quality scores
% N_obs - the number of ovservers
% C_mat - matrix with the comparisons to run in an experiment (incomplete
%         design). The size of the matrix must be length(q) x length(q).
%         The full design (all comparisons) are assumed if the matrix
%         is not provided. Only top-right half of the matrix should be filled -
%         the values will be copied to bottom-left half to ensure symmetry.
%         Example: to compare 1st and 2nd codition, put 1 in C_mat(1,2)
%         
%
% Author: Rafal Mantiuk


if( ~exist( 'pw_scale', 'file' ) )
    addpath( '../' );
end

% So that the result is the same each time
%s = RandStream('mt19937ar','Seed',0);
%RandStream.setGlobalStream(s);

N = length(q);

if( ~exist( 'C_mat', 'var' ) )
    C_mat = true(N,N);
else
    % Make copy elements from top-right to ensure the natrix is symmetrical
    for rr=1:N
        for cc=(rr+1):N
            C_mat(cc,rr) = C_mat(rr,cc);
        end
    end
end

sigma_cdf = 1.4826; % for this sigma normal cummulative distrib is 0.75 @ 1
sigma_q = sigma_cdf/sqrt(2);

N_reps = 1;

q_s = zeros(N_exps,N);

parfor ee=1:N_exps
    
    MM = zeros( N_obs, N*N );
    for oo=1:N_obs
        M = zeros(N,N);
        
        for rr=1:N_reps
            x_r = randn(1,N)*sigma_q + q;
            
            M = M + ((x_r' * ones(1,N)) > (ones(N,1) * x_r));
            
            M(~C_mat) = 0;
        end
        
        MM(oo,:) = M(:);
        
    end    
    q_s(ee,:) = pw_scale_bootstrp( MM, 1, { 'display', 'none' } );
    
end

end


