
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

all_conditions = unique( [D.condition_1] );

MM = zeros(length(OBSs), N*N);

n_scenes = length(SCs);

conditions_scene = zeros(n_scenes,numel(all_conditions));

% for each scene
for sc=1:n_scenes,
     
    display( sprintf( 'Scene: %s', SCs{sc} ) );
    
    Ds = D( strcmp( D.scene, SCs{sc} ), :);
    
    OBSs_scene = unique( Ds.observer );
    
    %conditions_scene{sc} = unique([Ds.condition_1]);
        
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
    end
    
end

% Any observer with a dist_L > 1.5 can be considered an outlier
[L,dist_L] = pw_outlier_analysis( MM );
    

