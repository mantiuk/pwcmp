function pw_metric_pairwise_correlation( S, M, group_column, condition_column, opt )
arguments
    S  % cell array with scaled subjective data, returned by pw_scale_table
    M
    group_column {mustBeText} = 'scene'
    condition_column {mustBeText} = 'condition'
    opt.alpha {mustBeInRange(opt.alpha,0,1)} = 0.05
    opt.reference_Q {mustBeNumeric} = 10;
end


% Collect all boostrap samples and scaled JODs into a single array


N_cond = numel(S{1}.conditions);

has_all = false;
for gg=1:length(S)
    if strcmp(S{gg}.group, 'all' )
        has_all = true;
        break;
    end
end
N_group = length(S) - double(has_all);
N_bstrp = size(S{1}.stats.bstrp,1);

bstrp = zeros(N_bstrp,N_cond*N_group);
jod = zeros(N_cond*N_group,1);
grps = cell(N_group,1);
conds = S{1}.conditions;

pp = 1;
for gg=1:length(S)
    group = S{gg}.group;
    if strcmp(group, 'all' )
        continue;
    end
    grps{pp} = group;
    pos = (pp-1)*N_cond+1;
    bstrp(:,pos:(pos+N_cond-1)) = S{gg}.stats.bstrp;
    jod(pos:(pos+N_cond-1)) = S{gg}.jod;
    pp = pp+1;
end

% Match metric quality score with the subjective dataset score
met_q = ones(N_cond*N_group,1)*opt.reference_Q; % Metric quality scores
for kk=1:height(M)
    ind_grp = find(strcmp(M.group{kk},grps));
    if isempty(ind_grp)
        error( 'Group "%s" present in the metric table, but missing in the subjective dataset', M.group{kk} );
    end
    ind_cond = find(strcmp(M.condition{kk},conds));
    if isempty(ind_grp)
        error( 'Condition "%s" present in the metric table, but missing in the subjective dataset', M.condition{kk} );
    end
    ind = (ind_grp-1)*N_cond + ind_cond;
    met_q(ind) = M.Q(kk);
end

% Find only significant pairs

jod_cov = cov(bstrp);
var_c = diag(jod_cov);
var_pair = var_c + var_c' - 2*jod_cov;

D = abs(jod-jod');

p = normcdf( -abs(D), 0, sqrt(var_pair) );
different = triu(p<(opt.alpha/2),1); % alpha/2 because we need a 2-tailed test
N = N_group*N_cond;
pair_total = N*(N-1)/2;
pair_diff = nnz(different);
fprintf( 1, 'Group: %s: %d out of %d pairs are statistically different (%.1f%%)\n', group, pair_diff, pair_total, pair_diff/pair_total*100 )

% for kk=1:(N-1)
%     Q_subj_D = vertcat( Q_subj_D, (Rs{gg}.jod(kk) - Rs{gg}.jod(different(kk,:))) );
%     Q_met_D = vertcat( Q_met_D, (Rs{gg}.Q(kk) - Rs{gg}.Q(different(kk,:))) );
%     grp = cell(nnz(different(kk,:)),1);
%     grp(:) = { Rs{gg}.group };
%     Q_group = vertcat( Q_group, grp );
% end




return

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