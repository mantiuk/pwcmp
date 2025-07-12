function [rho, rho_ci] = pw_metric_pairwise_correlation( S, M, within_set, opt )
arguments
    S  % cell array with scaled subjective data, returned by pw_scale_table
    M {istable(M)}
    within_set {mustBeText} = ''
    opt.group_column {mustBeText} = 'group'
    opt.condition_column {mustBeText} = 'condition'
    opt.alpha {mustBeInRange(opt.alpha,0,1)} = 0.05
    opt.scatter_plot {mustBeNumericOrLogical} = false
end
% The function computes the correlation between subjective scores (cell array
% S) and metric predictions (table M) for the pairs of conditions that are
% statistically different. The correlation is computed for the differences
% between the pairs (A,B) as corr( JOD_A - JOD_B, MET_A - MET_B ).
%
% This analysis offers two advantages:
% - we ignore the pairs for which there is insufficient statistical evidence
% to claim that one condition is better than another (the subjective data 
% is too noisy)
% - we can measure the correlation while discounting the effect of one or
% more factors. For example, we can check whether metric can predict the
% content-dependent differences in quality while ignoring all other factors. 
% This is controlled by the `within_set` parameter.
% 
% S - a cell array returned as a second return argument of
% pw_scale_table() with bootstrapping enabled. It is recommended that at
% least 500 bootstrap samples are used. It is good idea to store the cell
% array as a .mat file as running pw_scale_table() with bootstrapping can
% take some time.
%
% M - is a table with at least three columns:
% * 'group' - the name of the group used to group the conditions while scaling.
% This should correspond to the values in the column specified as group_col
% parameter of pw_scale_table(). Typically, it is a content_id.
% * 'condition' - the name of the conditions, which correspond to the values
% in the columns passed as `condition_cols` parameter of pw_scale_table().
% Typically, these are distortion types and distortion levels (combined as
% one string label).
% * 'Q'  - the quality prediction of a quality metric
%
% Note that table M must contain the same number of conditions and groups
% as in the subjective data (S), including reference conditions. 
%
% within_set - the name of the column of 'M', which is used to group pair
% of conditions into sets. The conditions are paired with each other within 
% each set, but not across the sets. For example, if you want to check 
% whether a metric can predict the effect of condition while ignoring the effect of
% the group (content), pass 'group'. This parameter shouyld be the name of 
% any column in 'M'.

if ~isempty(within_set) && ~ismember( within_set, M.Properties.VariableNames )
    error( 'Within group must be the column name in "M"')
end

% Collect all boostrap samples and scaled JODs into a single array
N_cond_all = 0;
N_group  = 0;
for gg=1:length(S)
    if strcmp(S{gg}.group, 'all' )
        continue;
    end
    N_cond_all = N_cond_all + numel(S{gg}.conditions);
    N_group = N_group+1;
end
N_bstrp = size(S{1}.stats.bstrp,1);

if height(M) ~= N_cond_all
    warning( 'The number of conditions in the table "M" (%d) must be the same as the number of conditions in the scaled subjective results "S" (%d). If table "M" is missing reference conditions, those must be added.', height(M), N_cond_all );
end

bstrp = zeros(N_bstrp,N_cond_all);
jod = zeros(N_cond_all,1);
grps = cell(N_group,1);
conds = S{1}.conditions;
grp_index = nan(N_group,1);

pp = 1;
pg = 1;
for gg=1:length(S)
    group = S{gg}.group;
    if strcmp(group, 'all' )
        continue;
    end
    grps{pg} = group;
    grp_index(pg) = pp;
    N_cond = numel(S{gg}.conditions);
    pos = pp;
    bstrp(:,pos:(pos+N_cond-1)) = S{gg}.stats.bstrp;
    jod(pos:(pos+N_cond-1)) = S{gg}.jod;
    pp = pp+N_cond;
    pg = pg+1;
end

% Match metric quality score with the subjective dataset score
met_q =  nan(N_cond_all,1); % Metric quality scores
M.index = nan(height(M),1);
for kk=1:height(M)
    ind_grp = find(strcmp(M.(opt.group_column){kk},grps));
    if isempty(ind_grp)
        error( 'Group "%s" present in the metric table "M", but missing in the scaled subjective data "S".', M.group{kk} );
    end
    ind_cond = find(strcmp(M.(opt.condition_column){kk},S{ind_grp}.conditions));
    if isempty(ind_cond)
        error( 'Condition "%s" present in the metric table "M", but missing in the scaled subjective data "S".', M.(opt.condition_column){kk} );
    end
    ind = grp_index(ind_grp) + ind_cond-1;
    met_q(ind) = M.Q(kk);    
    M.index(kk) = ind;
end
assert( all(~isnan(met_q)) );
assert( all(~isnan(M.index)) );

if isempty( within_set )
    cor_sets = { 'all' };
else
    % consider pairs differences only within the groups
    cor_sets = unique( M.(within_set) );
end

Q_set = [];
pairs = [];

for gg=1:length(cor_sets)

    if isempty( within_set )
        inds = 1:size(bstrp,2);
    else
        inds = M.index(strcmp(M.(within_set),cor_sets{gg}));
    end

    % Find only significant pairs
    jod_cov = cov(bstrp(:,inds));
    var_c = diag(jod_cov);
    var_pair = var_c + var_c' - 2*jod_cov;
    
    D = abs(jod(inds)-jod(inds)');
    
    p = normcdf( -abs(D), 0, sqrt(var_pair) );
    different = triu(p<(opt.alpha/2),1); % alpha/2 because we need a 2-tailed test
    N = numel(inds);
    pair_total = N*(N-1)/2;
    pair_diff = nnz(different);
    fprintf( 1, 'Set: %s: %d out of %d pairs are statistically different (%.1f%%)\n', cor_sets{gg}, pair_diff, pair_total, pair_diff/pair_total*100 )

    new_pairs = zeros(pair_diff,2);

    pp_beg = 1;
    for kk=1:(N-1)
        pp_end = pp_beg + nnz(different(kk,:)) -1;
        new_pairs(pp_beg:pp_end,1) = inds(kk);
        new_pairs(pp_beg:pp_end,2) = inds(different(kk,:));
        pp_beg = pp_end+1;
    end
    q_set = cell(pair_diff,1);
    q_set(:) = { cor_sets{gg} };
    Q_set = vertcat( Q_set, q_set );
    pairs = cat( 1, pairs, new_pairs );
    
end

Q_subj_D = jod(pairs(:,1)) - jod(pairs(:,2));
Q_met_D = met_q(pairs(:,1)) - met_q(pairs(:,2));

rho = corr( Q_met_D, Q_subj_D, 'Type', 'Pearson' );

% Use boostrapping samples to estimate the confidence interval of the
% correlations
if true
    % Random sample
    N_samples = N_bstrp*floor(sqrt(N_bstrp));
    brs = randi(N_bstrp, N_samples, 2);
else
    % Every combination of boostrap samples - slow, but most accurate
    [b1, b2] = ind2sub( [N_bstrp, N_bstrp], 1:N_bstrp^2 );
    brs = cat( 2, b1', b2' );
end

% Use boostrap samples
B_subj_D = bstrp(brs(:,1),pairs(:,1)) - bstrp(brs(:,2),pairs(:,2));

rho_dist = corr( Q_met_D, B_subj_D', 'Type', 'Pearson' );
rho_ci = prctile( rho_dist', [opt.alpha/2*100, 100-opt.alpha/2*100], 'all' );

if opt.scatter_plot
    clf;
    plot( [0 0], [min(Q_subj_D) max(Q_subj_D)], '--k' );
    hold on
    plot( [min(Q_met_D) max(Q_met_D)], [0 0], '--k' );
    if length(Q_set)>1 && length(Q_set)<20
        gscatter( Q_met_D, Q_subj_D, Q_set );
    else
        scatter( Q_met_D, Q_subj_D, 'o', MarkerEdgeColor='b', MarkerFaceColor='b' );
    end
    title( sprintf( 'PLCC=%g (%g, %g)', rho, rho_ci(1), rho_ci(2) ) );
    xlabel( 'Metric A-B' )
    ylabel( 'Subjective score A-B [JOD]' )
end


end