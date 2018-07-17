function [bias, ci, d, q_m, meas_mse] = simulate_exp_bootstrp( q, N_obs, do_ties, full_design )

if( ~exist( 'pw_scale', 'file' ) )
    addpath( '../' );
end

% So that the result is the same each time
%s = RandStream('mt19937ar','Seed',0);
%RandStream.setGlobalStream(s);

N = length(q);

sigma_cdf = 1.4826; % for this sigma normal cummulative distrib is 0.75 @ 1
sigma_q = sigma_cdf/sqrt(2);

N_exps = 10000;
N_reps = 1;

q_s = zeros(N_exps,N);

for ee=1:N_exps
    
    MM = zeros( N_obs, N*N );
    for oo=1:N_obs
        M = zeros(N,N);
        
        for rr=1:N_reps
            x_r = randn(1,N)*sigma_q + q;
            
            if( do_ties )
                tie_thr = max( 0, 0.7+randn()*0.3 );
                D = ((x_r' * ones(1,N)) - (ones(N,1) * x_r));
                M = M + (D>=tie_thr);
                M = M + ((D>-tie_thr & D<tie_thr)*0.5).*~eye(size(D));
            else
                M = M + ((x_r' * ones(1,N)) > (ones(N,1) * x_r));
            end
            if( ~full_design )
                % Only neighbours
                ss = (triu( ones(N), 2 ) | tril( ones(N), -2 ));
                M(ss) = 0;
            end
        end
        
        MM(oo,:) = M(:);
        
    end
    
    q_s(ee,:) = pw_scale_bootstrp( MM, 1, { 'display', 'none' } );
    
end

q_m = mean(q_s);
e_up = prctile(q_s,97.5)-q_m;
e_low = q_m - prctile(q_s,2.5);

res_true_qm = q_m-q;
ci = mean( cat( 2, e_up(2:end), e_low(2:end) ) );
bias = mean( res_true_qm(2:end) );
res_est = q_s-repmat( q_m, [size(q_s,1) 1] );

d = mean( diff(q_m)./std( res_est(:,2:end) ) );
meas_mse = mean( res_true_qm(2:end).^2 );

end

