function [L,dist_L] = pw_outlier_analysis_table( T, group_col, condition_cols, observer_col, selection_col )
arguments
    T table
    group_col char
    condition_cols (1,2) cell = { 'condition_A', 'condition_B' }
    observer_col char = 'observer'
    selection_col char = 'selected_A'
end
% Performs outlier analysis on the data in a table. The observers for which
% the log-likelihood is significantly smaller are potential outliers. 
%
% See pw_outlier_analysis for more information on the outlier analysis. The
% method is explained in Section 9 of https://arxiv.org/abs/1712.03686
% 
% [L,dist_L] = pw_outlier_analysis_table( T, group_col, condition_cols, observer_col, selection_col )
%
% T is the table with the results
% group_col - name of the column used to group the results. For example, if
%             the scaling is performed separately for each scene/content, 
%             this should be the column with the ID of the content.
% conditions_col - a cell array with the name of the column string two
%             compared conditions, for example { 'condition_A', 'condition_B' }
% observer_col - the name of the column storing the ids of the observers
% selection_col - the name of the variable indicating whether the first
%             condition was selected. The values must be 0 or 1

if isempty(group_col)
    GRs = {};
else
    GRs = unique( T.(group_col) ); % list of groups
end
C = unique( cat( 1, T.(condition_cols{1}), T.(condition_cols{2}) ) ); % all conditions

% The observer column must contain strings
if isnumeric(T.(observer_col))
    T.(observer_col) = mat2cell(num2str(T.(observer_col)), ones(height(T),1));
end

N = length(C);

%MM_l = {};

OBSs_all = unique( T.(observer_col) );

N_obs = numel(OBSs_all);
L_all = zeros(N_obs,length(GRs));

% for each group
for gg=1:length(GRs)

    Ds = T( strcmp( T.(group_col), GRs{gg} ), :);
    group = GRs{gg};
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
    
    [L,dist_L] = pw_outlier_analysis(MM);
    
    % Each group may have a different number of observers, so we need to
    % match those
    for oo=1:length( OBSs )
        ind = find( strcmp(OBSs{oo}, OBSs_all) );
        assert( numel(ind)==1 );
        L_all(ind,gg) = L(oo);
    end 
end


clf;
ind = 1:N_obs;
plot( sum(L_all,2), ind, 'o' );
set( gca, 'YTick', ind );
set( gca, 'YTickLabel', OBSs_all );
xlabel( 'Log likelihood' );
title( 'The log-likelihood reports similarity to all other observers', 'FontWeight', 'normal' );
grid on;

end