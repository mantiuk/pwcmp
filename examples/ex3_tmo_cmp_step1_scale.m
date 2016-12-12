% Part 1 of the ex3. The example shows how to process actual data from the
% experiment. The data is loaded from the comma-separated file
% 'ex3_tmo_cmp_data2.csv'. Then, the results are scaled in JOD units and
% the result is saved in two formats: in another comma separated file and
% in a .mat file.
%
% Since computing bootstrapping takes time, it is advisable to separate
% scaling from plotting. 

if( ~exist( 'pw_scale', 'file' ) )
    addpath( '../' );
end

bootstrap_samples = 500;

D = dataset( 'File', 'ex3_tmo_cmp_data2.csv', 'Delimiter', ',' );
D = D( ~strcmp( D.scene, 'peniches' ), : );

SCs = unique( D.scene ); % list of scenes
C = unique( cat( 1, D.condition_1, D.condition_2 ) ); % all conditions

N = length(C);

R = []; % Store scaled results as a dataset (for the CSV file)
Rs = cell(length(SCs),1); % and as a cell array

% for each scene
for sc=1:length(SCs)
    
    display( sprintf( 'Scene: %s', SCs{sc} ) );
    
    Ds = D( strcmp( D.scene, SCs{sc} ), :);
    
    OBSs = unique( Ds.observer );
    
    jod = zeros( N, 1 );
    
    boostrap_samples = 10;
    
    MM = zeros(length(OBSs), N*N);
    
    for oo=1:length( OBSs )
        
        Dso = Ds(strcmp( Ds.observer, OBSs{oo} ), :);
        
        M = zeros(N,N);
        
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
        
        MM(oo,:) = M(:);
        
    end
    
    tic
    [jod, stats] = pw_scale_bootstrp( MM, bootstrap_samples, { 'use_parallel', 'never' } );
    toc
    
    
    Rn = dataset( ...
        { repmat( { SCs{sc} }, [N 1] ), 'scene' }, { C, 'tmo' }, ...
        { jod, 'jod' }, { stats.jod_low, 'jod_low' }, { stats.jod_high, 'jod_high' } ...
        );
        
    if( ~isempty(R) )
        R = vertcat( R, Rn );
    else
        R = Rn;
    end    
    
    Rs{sc}.scene = SCs{sc};
    Rs{sc}.jod = jod;
    Rs{sc}.stats = stats;
    
end


export( R, 'file', 'ex3_tmo_cmp_scaled.csv', 'Delimiter', ',' );
save( 'ex3_tmo_cmp_scaled.mat', 'Rs', 'C' );
