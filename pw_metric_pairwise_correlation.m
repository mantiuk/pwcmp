function pw_metric_pairwise_correlation( S, M, group_column, condition_column, opt )
arguments
    S 
    M
    group_column {mustBeText} = 'scene'
    condition_column {mustBeText} = 'condition'
    opt.alpha {mustBeInRange(opt.alpha,0,1)} = 0.05
    opt.reference_Q {mustBeNumeric} = 10;
end

Rs = S;

Q_subj_D = [];
Q_met_D = [];
Q_group = [];

for gg=1:length(Rs)

    group = Rs{gg}.group;
    if strcmp(group, 'all' )
        continue;
    end
    Ms = M(strcmp(M.(group_column), group),:);
    N = length(Rs{gg}.jod);
    Rs{gg}.Q = ones(N,1) * opt.reference_Q;
    for kk=1:N
        ind = find(strcmp(Ms.(condition_column),Rs{gg}.conditions{kk}),1);
        if ~isempty(ind)
            Rs{gg}.Q(kk) = Ms.Q(ind);
        end
    end

    % Find only significant pairs
    var_c = diag(Rs{gg}.stats.jod_cov);
    var_pair = var_c + var_c' - 2*Rs{gg}.stats.jod_cov;

    D = abs(Rs{gg}.jod-Rs{gg}.jod');

    p = normcdf( -abs(D), 0, sqrt(var_pair) );
    different = triu(p<(opt.alpha/2),1); % alpha/2 because we need a 2-tailed test
    pair_total = N*(N-1)/2;
    pair_diff = nnz(different);
    fprintf( 1, 'Group: %s: %d out of %d pairs are statistically different (%.1f%%)\n', group, pair_diff, pair_total, pair_diff/pair_total*100 )
    
    for kk=1:(N-1)
        Q_subj_D = vertcat( Q_subj_D, (Rs{gg}.jod(kk) - Rs{gg}.jod(different(kk,:))) );
        Q_met_D = vertcat( Q_met_D, (Rs{gg}.Q(kk) - Rs{gg}.Q(different(kk,:))) );
        grp = cell(nnz(different(kk,:)),1);
        grp(:) = { Rs{gg}.group };
        Q_group = vertcat( Q_group, grp );
    end

end

figure(1);
clf;

gscatter( Q_met_D, Q_subj_D, Q_group );
rho = corr( Q_met_D, Q_subj_D, 'Type', 'Pearson' );
title( sprintf( 'PLCC=%g', rho ) );
xlabel( 'Metric A-B' )
ylabel( 'Subjective score A-B [JOD]' )

end