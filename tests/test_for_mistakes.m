% Simulating people making mistakes - Monte Carlo

if( ~exist( 'pw_scale', 'file' ) )
    addpath( '../' );
end

x = [0 1 2 3 4];  % Correct response for a typical observer

N_samples = 1000;

N_obs = 8;
N_rep = 3;

N_c = length(x);

p_mist = 0.05;

sigma_cdf = 1.4826; % for this sigma normal cummulative distrib is 0.75 @ 1

N = N_obs * N_rep * (N_c*(N_c-1)/2);

QQ = zeros( N_samples, N_c );
for ss=1:N_samples
    
    if mod( ss, 100 ) == 1
        fprintf( 1, '.' );
    end
    MM = zeros(N_obs,N_c*N_c);
    for oo=1:N_obs  % For each observer
        M = zeros(N_c);
        for rr=1:N_rep
            for ii=1:(N_c-1)
                for jj=(ii+1):N_c
                    
                    rm = rand(1);
                    if rm <= p_mist % Random mistake
                        if( rm <= (p_mist/2) )
                            M(jj,ii) = M(jj,ii) + 1;
                        else
                            M(ii,jj) = M(ii,jj) + 1;
                        end
                    else
                        d = x(jj)-x(ii);
                        rp = randn()*sigma_cdf;
                        if( rp <= d )
                            M(jj,ii) = M(jj,ii)+1;
                        else
                            M(ii,jj) = M(ii,jj) + 1;
                        end
                    end
                end
            end
        end
        MM(oo,:) = M(:);
    end
    
    
    D = reshape( sum(MM), [N_c N_c] );
    
    Q = pw_scale( D );
    
    QQ(ss,:) = Q';
    
end

%%

clf
plot( x, x, '--k' );
hold on
errorbar( x, mean(QQ), std(QQ), 'x' );
hold off

xlim( [min(x)-0.1 max(x)+0.1] );
xlabel( 'True quality scores' );
ylabel( 'Scaling results' );

MSE =mean( mean( (QQ - repmat( x, [size(QQ,1) 1] )).^2 ) );
rep_text = sprintf( 'bias=%g  MSE=%g', mean( x-mean(QQ) ), MSE );
title( rep_text );
