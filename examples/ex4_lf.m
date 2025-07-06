% This is an example of scaling a large-size experiment (26580 comparisons) 
%
% Running this example can take some time (up to a few hours). Set
% boostrap_samples = 0 to compute JODs but no confidence intervals. 
%
% The data for this example comes from the paper:
% Adhikarla, V. K., Vinkler, M., Sumin, D., Mantiuk, R. K., Myszkowski, K., Seidel, H.-P., & Didyk, P. (2017). 
% Towards a quality metric for dense light fields. 
% Computer Vision and Pattern Recognition (CVPR), 58–67.
%
% Project web page: http://lightfields.mpi-inf.mpg.de/

if( ~exist( 'pw_scale', 'file' ) )
    addpath( '../' );
end

D = dataset( 'File', 'lf_quality_pwcmps.csv', 'Delimiter', ',' );

STs = unique( D.session );
SCs = unique( D.scene );

% A condition is a combination of distotion type and distortion level. 
% Create a column for each compared condition.
D.condition_1 = strcat( D.dist_type1, '_', num2str( D.dist_level1 ) );
D.condition_2 = strcat( D.dist_type2, '_', num2str( D.dist_level2 ) );

R = [];

boostrap_samples = 500;

% for each session type
for st=1:length(STs)
    
    % for each scene
    for sc=1:length(SCs)
        
        fprintf( 1, 'Scene: %s\n', SCs{sc} );
        
        Ds = D(strcmp( D.session, STs{st} ) & strcmp( D.scene, SCs{sc} ), :);

        C = unique( cat( 1, Ds.condition_1, Ds.condition_2 ) ); % all conditions (for this scene)
        N = length(C);
        
        OBSs = unique( Ds.observer );
        
        %jnd = zeros( N, 1 );
                
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
        [jod, stats] = pw_scale_bootstrp( MM, boostrap_samples , 'regularization', 'fix0', 'prior', 'gaussian', 'use_parallel', 'always' );
        toc
        jod_se = sqrt(diag(stats.jod_cov));
        
        dist_level = zeros( N, 1 );
        dist_type = cell( N, 1 );
        for kk=1:N
            if( strcmp( C{kk}, 'original' ) )
                dist_type{kk} = 'original';
                dist_level(kk) = 0;
            else
                F = textscan( C{kk}, '%s%d', 'Delimiter', '_' );
                dist_type{kk} = F{1}{1};
                dist_level(kk) = F{2};
            end
        end
        
        Rn = dataset( { dist_type, 'dist_type' }, { dist_level, 'dist_level' }, ...
            { jod, 'jod' }, { stats.jod_low, 'jod_low' }, { stats.jod_high, 'jod_high' }, ...
            { jod_se, 'jod_se' }, ...
            { repmat( { STs{st} }, [N 1] ), 'session_type' }, ...
            { repmat( { SCs{sc} }, [N 1] ), 'scene' } );
                
        if( ~isempty(R) )
            R = vertcat( R, Rn );
        else
            R = Rn;
        end
        
    end
    
end

export( R, 'file', 'lf_quality_scaled.csv', 'Delimiter', ',' );

