
if( ~exist( 'pw_scale', 'file' ) )
    addpath( '../' );
end

bootstrap_samples = 500;

D = dataset( 'File', '/home/mp867/Desktop/pairwise_cmp/2017_compression_colorspace/colorCompressionExperiment_subjectivePCdata_comma.csv', 'Delimiter', ',' );
D = D( ~strcmp( D.scene, 'peniches' ), : );

SCs = unique( D.scene ); % list of scenes
C = unique( cat( 1, D.condition_1, D.condition_2 ) ); % all conditions

N = length(C);

R = []; % Store scaled results as a dataset (for the CSV file)
Rs = cell(length(SCs),1); % and as a cell array
all = zeros(N,N);
    
OBSs = unique( D.observer );

MM = zeros(length(OBSs), N*N);

n_scenes = length(SCs);

c_scene = zeros(n_scenes,numel(C));

% We need to define first which is the outlier in our dataset
outlier = 7;

% for each scene
for sc=1:n_scenes,
     
    display( sprintf( 'Scene: %s', SCs{sc} ) );
    
    Ds = D( strcmp( D.scene, SCs{sc} ), :);
    
    OBSs_scene = unique( Ds.observer );
    
    conditions = unique([Ds.condition_1 Ds. condition_2]);
    
    for cc=1:length( conditions ),
        
       c_scene(sc,:) = c_scene(sc,:) + (strcmp( C, conditions{cc} )');
        
    end
        
    for oo=1:length( OBSs_scene )
        
        Dso = Ds(strcmp( Ds.observer, OBSs_scene{oo} ), :);
        
        M = zeros(N,N);
                    
        if min(Dso.selection)==1,
            Dso.selection = Dso.selection - 1;
        end
        
        for kk=1:length(Dso)
            
            % Find the indexes of both conditions
            c1 = find( strcmp( Dso.condition_1(kk), C ), 1 );
            c2 = find( strcmp( Dso.condition_2(kk), C ), 1 );
            
            if( isempty( c1 ) )
                error( 'Cannot find condition %s', Dso.condition_1(kk) );
            end
            if( isempty( c2 ) )
                error( 'Cannot find condition %s', Dso.condition_2(kk) );
            end
            
            if( Dso.selection(kk) == 0 )
                M(c1,c2) = M(c1,c2) + 1; % c1 better than c2
            else
                M(c2,c1) = M(c2,c1) + 1; % c2 better than c1
            end
        end
        
        % find that observer (different observers can do different
        % experiments)
        index = find(strcmp(OBSs, OBSs_scene{oo})==1);
        MM(index,:) = MM(index,:) + M(:)';
        all = all + M;
            
        prob = M./(M + M');
        prob(isnan(prob)) = 0;
        
        s_columns = sum(prob,1);
        s_rows = sum(prob,2);
        
        data_boxplot_columns(oo,:) = (s_columns./(s_columns + s_rows'));
        data_boxplot_rows(oo,:) = s_rows./(s_columns' + s_rows);        
    end
    
% Save the data from the outlier
o_columns = data_boxplot_columns(outlier,:);
o_rows = data_boxplot_rows(outlier,:);    

% Remove the data from the outlier from the distribution
data_boxplot_columns(outlier,:) = [];
data_boxplot_rows(outlier,:) = [];

subplot(3,2,sc);
hold on;
% Plot the data for that scene
boxplot(data_boxplot_rows(:,c_scene(sc,:)==1),C(c_scene(sc,:)==1))
plot([1:(sum(c_scene(sc,:)))],o_rows(c_scene(sc,:)==1),'ok','Linewidth',2)
xlabel('Conditions compared')
ylabel('Probability of higher quality')
legend('Outlier')

title(SCs{sc})
hold off

end
