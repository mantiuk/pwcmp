% This script tests pw_outlier_analysis function
% It simulates an experiment in which one observer (observer 1) gives
% responses that are different than the responses of all other observers

if( ~exist( 'pw_scale', 'file' ) )
    addpath( '../' );
end

x = [0 1 2 3 4];  % Correct response for a typical observer
x_outl = [0 1 2 3 4] * -1;  % False response for an outlier

N_obs = 20;
N_rep = 2;

N_scenes = 5;

N_c = length(x);

sigma_cdf = 1.4826; % for this sigma normal cummulative distrib is 0.75 @ 1

N = N_obs * N_rep * (N_c*(N_c-1)/2);

M_list = cell(N_scenes,1);
for ss=1:N_scenes % For each scene
    
    MM = zeros(N_obs,N_c*N_c);
    for oo=1:N_obs  % For each observer
        M = zeros(N_c);
        for rr=1:N_rep
            for ii=1:(N_c-1)
                for jj=(ii+1):N_c
                    
                    if oo == 1
                        d = x_outl(jj)-x_outl(ii);
                    else
                        d = x(jj)-x(ii);
                    end
                    rp = randn()*sigma_cdf;
                    if( rp <= d )
                        M(jj,ii) = M(jj,ii)+1;
                    else
                        M(ii,jj) = M(ii,jj) + 1;
                    end
                end
            end
        end
        MM(oo,:) = M(:);
    end
    
    M_list{ss} = MM;
end

L = pw_outlier_analysis( M_list );

clf
bar( -L );
xx = 1:length(L);
set( gca, 'XTick', xx );
set( gca, 'XTickLabel', num2str( xx' ) );

