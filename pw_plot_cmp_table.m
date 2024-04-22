function pw_plot_cmp_table( T, group_col, condition_cols, observer_col, selection_col )
arguments
    T table
    group_col char
    condition_cols (1,2) cell = { 'condition_A', 'condition_B' }
    observer_col char = 'observer'
    selection_col char = 'selected_A'
end
% Plot a table with pairwise comparisons. 
%
% pw_plot_cmp_table( T, group_col, condition_cols, observer_col, selection_col )
%
% T is the table with the results
% group_col - name of the column used to group the results. For example, if
%             the scaling is performed separately for each scene/content, 
%             this should be the column with the ID of the content. Pass
%             an empty array ([]) if there are no groups in the table. 
% conditions_col - a cell array with the name of the column string two
%             compared conditions, for example { 'condition_A', 'condition_B' }
% observer_col - the name of the column storing the ids of the observers
% selection_col - the name of the variable indicating whether the first
%             condition was selected. The values must be 0 or 1

if isempty(group_col)
    GRs = {};
    N_g = 1;
else
    GRs = unique( T.(group_col) ); % list of groups
    N_g = length( GRs );
end
C = unique( cat( 1, T.(condition_cols{1}), T.(condition_cols{2}) ) ); % all conditions

% The observer column must contain strings
if isnumeric(T.(observer_col))
    T.(observer_col) = mat2cell(num2str(T.(observer_col)), ones(height(T),1));
end

N = length(C);


% for each group
for gg=1:N_g

    if isempty(group_col)
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
    
    cmp_table = reshape( sum(MM), N, N );

    figure(gg);
    clf;
    C_noun = strrep( C, '_', ' ' );
    heatmap( C_noun, C_noun, cmp_table );
    title( 'The condition in a row selected over the condition in a column' );
 
end

end