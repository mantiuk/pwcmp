function [R, Rs] = pw_scale_table( T, group_col, condition_cols, observer_col, selection_col, options )
arguments
    T table
    group_col char
    condition_cols (1,2) cell = { 'condition_A', 'condition_B' }
    observer_col char = 'observer'
    selection_col char = 'selected_A'
    options.bootstrap_samples (1,1) {mustBeGreaterThanOrEqual(options.bootstrap_samples,0),mustBeInteger} = 500
    options.prior char = 'gaussian'
    options.regularization char {mustBeMember(options.regularization, {'mean0', 'fix0'})} = 'mean0'
    options.do_all logical = false
    options.anchor_condition {mustBeText} = ''
end
% Scales pairwise comparison results stored in a table
% 
% [R, Rs] = pw_scale_table( T, group_col, condition_cols, observer_col, selection_col )
% [R, Rs] = pw_scale_table( T, ..., 'boostrap_samples', N )
% [R, Rs] = pw_scale_table( T, ..., 'prior', prior_name )
% [R, Rs] = pw_scale_table( T, ..., 'do_all', bool )
%
% T is the table with the results
% group_col - name of the column used to group the results. The scaling is performed
%             separately on each group. Pass an empty string if all the
%             data should be scaled together without splitting into
%             groups.
% conditions_col - a cell array with the name of the column string two
%             compared conditions, for example { 'condition_A', 'condition_B' }
% observer_col - the name of the column storing the ids of the observers
% selection_col - the name of the variable indicating whether the first
%             condition was selected. The values must be 0 or 1
% optional arguments:
% bootstrap_samples - see 'bootstrap_samples' in pw_scale_bootstrp
% prior - see 'prior' in pw_scale_bootstrp
% do_all - in addition to the per-group scaling, scale all results across
%        all the groups

if isempty(group_col)
    GRs = {};
    options.do_all = true;
else
    GRs = unique( T.(group_col) ); % list of groups
end
C = unique( cat( 1, T.(condition_cols{1}), T.(condition_cols{2}) ) ); % all conditions

if ~isempty( options.anchor_condition )
    ind = find( strcmp( C, options.anchor_condition ) );
    if length(ind) ~= 1
        error( 'Anchor condition %s not found', options.anchor_condition );
    end
    C{ind} = C{1};
    C{1} = options.anchor_condition;
end

% The observer column must contain strings
if isnumeric(T.(observer_col))
    T.(observer_col) = mat2cell(num2str(T.(observer_col)), ones(height(T),1));
end

N = length(C);

if options.do_all
    start_group = 0;
else
    start_group = 1;
end

R = []; % Store scaled results as a dataset (for the CSV file)
Rs = cell(length(GRs) + 1-start_group,1); % and as a cell array

pp = 1;
% for each scene
for gg=start_group:length(GRs)

    if gg==0
        Ds = T;
        group = 'all';
    else
        Ds = T( strcmp( T.(group_col), GRs{gg} ), :);
        group = GRs{gg};
    end
    fprintf( 1, 'Group: %s\n', group );
    
    OBSs = unique( Ds.(observer_col) );
        
    MM = zeros(length(OBSs), N*N);
    
    for oo=1:length( OBSs )
        
        Dso = Ds(strcmp( Ds.(observer_col), OBSs{oo} ), :);
        
        M = zeros(N,N);
        
        for kk=1:height(Dso)
            
            % Find the indexes of both conditions
            c1 = find( strcmp( Dso.(condition_cols{1})(kk), C ), 1 );
            c2 = find( strcmp( Dso.(condition_cols{2})(kk), C ), 1 );
            
            if( isempty( c1 ) )
                error( 'Cannot find condition %s', Dso.(condition_cols{1})(kk) );
            end
            if( isempty( c2 ) )
                error( 'Cannot find condition %s', Dso.(condition_cols{2})(kk) );
            end
            
            if( Dso.(selection_col)(kk) )
                M(c1,c2) = M(c1,c2) + 1; % c1 better than c2
            else
                M(c2,c1) = M(c2,c1) + 1; % c2 better than c1
            end
        end
        
        MM(oo,:) = M(:);
        
    end
    
    tic
    [jod, stats] = pw_scale_bootstrp( MM, options.bootstrap_samples, 'prior', options.prior, 'regularization', options.regularization );
    toc    
    
    if ~isempty(group_col)
        S.(group_col) = repmat( { group }, [N 1] );
    end

    if ~isempty( options.anchor_condition ) % 
        offset = jod(1); % Shift the JODs so that the anchor=0
    else
        offset = 0;
    end

    S.condition = C;
    S.jod = jod - offset;
    S.jod_low = stats.jod_low - offset;
    S.jod_high = stats.jod_high - offset;

    Rn = struct2table(S);
        
    if( ~isempty(R) )
        R = vertcat( R, Rn );
    else
        R = Rn;
    end    
    
    Rs{pp}.conditions = C;
    Rs{pp}.group = group;
    Rs{pp}.jod = jod;
    Rs{pp}.stats = stats;
    pp = pp + 1;
    
end

end