if( ~exist( 'pw_scale', 'file' ) )
    addpath( '../' );
end

D = dataset( 'File', 'db_subject_pw.csv', 'Delimiter', ',' );
SCs = unique( D.scene );

% Mark all Skip0 as reference
for kk=1:length(D)
    if( ~isempty( strfind( D.condition_1{kk}, 'Skip0' ) ) )
        D.condition_1{kk} = '0ref';
    end
    if( ~isempty( strfind( D.condition_2{kk}, 'Skip0' ) ) )
        D.condition_2{kk} = '0ref';
    end
end

C = unique( cat( 1, D.condition_1, D.condition_2 ) ); % all conditions

N = length(C);

R = [];

PC = struct();
PC.C = C;

% for each scene
for sc=1:length(SCs)
    
    display( sprintf( 'Scene: %s', SCs{sc} ) );
    
    Ds = D( strcmp( D.scene, SCs{sc} ), :);
    
    OBSs = unique( Ds.observer );
    
    jnd = zeros( N, 1 );
    
    boostrap_samples = 10;
    
    MM = zeros(length(OBSs), N*N);
    
    for oo=1:length( OBSs )
        
        Dso = Ds(strcmp( Ds.observer, OBSs{oo} ), :);
        
        M = zeros(N,N);
        
        for kk=1:length(Dso)
            c1 = find( strcmp( Dso.condition_1(kk), C ), 1 );
            c2 = find( strcmp( Dso.condition_2(kk), C ), 1 );
            
            if( isempty( c1 ) || isempty( c2 ) )
                warning( 'Cannot find condition' );
                continue;
            end
            
            if( Dso.selected(kk) == 1 )
                M(c1,c2) = M(c1,c2) + 1;
            else
                M(c2,c1) = M(c2,c1) + 1;
            end
        end
        
        MM(oo,:) = M(:);
        
    end
    
    tic
    [jnd, stats] = pw_scale_bootstrp( MM, 50 );
    toc
    jnd_se = sqrt(diag(stats.jnd_cov));
    
    dist_level = zeros( N, 1 );
    dist_type = cell( N, 1 );
    for kk=1:N
        if( strcmp( C{kk}, '0ref' ) )
            dist_type{kk} = 'ref';
            dist_level(kk) = 0;
        else
            F = textscan( C{kk}, '%sSkip%d', 'Delimiter', '_' );
            dist_type{kk} = F{1}{1};
            dist_level(kk) = F{2};
        end
    end
    
    Rn = dataset( { dist_type, 'dist_type' }, { dist_level, 'dist_level' }, ...
        { jnd, 'jnd' }, { stats.jnd_low, 'jnd_low' }, { stats.jnd_high, 'jnd_high' }, ...
        { jnd_se, 'jnd_se' }, ...
        { repmat( { SCs{sc} }, [N 1] ), 'scene' } );
    
    PC.(SCs{sc}) = stats.jnd_cov;
    
    if( ~isempty(R) )
        R = vertcat( R, Rn );
    else
        R = Rn;
    end
    
end


%export( R, 'file', 'pw_scaled.csv', 'Delimiter', ',' );
%save( 'pw_covariance.mat', 'PC' );

