function simulate_exp( x, N_obs, N_rep )
% Simulate experimemnt and run scaling
%
% x - true JND values for the conditions
% N_obs - number of observers
% N_rep - number of repetitions


% So that the result is the same each time
%s = RandStream('mt19937ar','Seed',0);
%RandStream.setGlobalStream(s);

N_c = length(x);

sigma_cdf = 1.4826; % for this sigma normal cummulative distrib is 0.75 @ 1

%N_obs = 30;
%N_rep = 4;

N = N_obs * N_rep * (N_c*(N_c-1)/2);

MM = zeros(N_obs,N_c*N_c);
for oo=1:N_obs
    M = zeros(N_c);
    for rr=1:N_rep
        for ii=1:(N_c-1)
            for jj=(ii+1):N_c
                                
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
    MM(oo,:) = M(:);
end

%M_total = reshape( sum(MM), size(M) );
%P = M_total ./ (M_total+M_total');
%D = norminv( P, 0, sigma_cdf );

[jnd, stats] = pw_scale_bootstrp( MM, 500 );

clf
html_change_figure_print_size( gcf, 8, 8 );
le = jnd-stats.jnd_low;
ue = stats.jnd_high-jnd;

max_xx = max(x)+2;
xx = linspace( 0, max_xx );
figure(1);
plot( xx, xx, '--k' );
hold on
errorbar( x, jnd, le, ue, 'o' );
hold off
grid on

xlim( [0 max_xx] );
ylim( [0 max_xx] );
xlabel( 'True mean quality score' );
ylabel( 'Measured quality score' );
set( gca, 'XTick', 0:max_xx );
set( gca, 'YTick', 0:max_xx );

figure(2);
N = length(jnd);
CONDs = cell(N,1);
for kk=1:N
    CONDs{kk} = num2str( x(kk) );
end
pw_plot_ranking_triangles( jnd, stats, CONDs );


end
