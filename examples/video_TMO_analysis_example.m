% Analyse the data in the video_TMO example: outlier detection, scaling and
% statistical analysis

if( ~exist( 'pw_scale', 'file' ) )
    addpath( genpath('../') );
end

bootstrap_samples = 500;

D = dataset( 'File', './ex3_tmo_cmp_data.csv', 'Delimiter', ',' );
D = D( ~strcmp( D.scene, 'peniches' ), : );

SCs = unique( D.scene ); % list of scenes
C = unique( cat( 1, D.condition_1, D.condition_2 ) ); % all conditions
% set TMO as first condition
C = circshift(C,1);

N = length(C);

R = []; % Store scaled results as a dataset (for the CSV file)
Rs = cell(length(SCs),1); % and as a cell array
all = zeros(N,N);
    
OBSs = unique( D.observer );

MM = zeros(length(OBSs), N*N);

n_scenes = length(SCs);

% for each scene
for sc=1:n_scenes
     
    fprintf( 1, 'Scene: %s\n', SCs{sc} );
    
    Ds = D( strcmp( D.scene, SCs{sc} ), :);
    
    OBSs_scene = unique( Ds.observer );
            
    % for each observer in that scene
    for oo=1:length( OBSs_scene )
        
        Dso = Ds(strcmp( Ds.observer, OBSs_scene{oo} ), :);
        
        M = zeros(N,N);
                    
        if min(Dso.selection)==1
            Dso.selection = Dso.selection - 1;
        end
        
        % read comparisons made by observer
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

% Simulate outlier observer - answers completely at random!
%aux = triu((randi(2,[N,N])-1),1);
%rand_obs = aux + triu(4-aux,1)';
%MM(index+1,:) = rand_obs(:);

% Any observer with a dist_L > 1.5 can be considered an outlier
[L,dist_L] = pw_outlier_analysis( MM );

% Find outliers
outliers = find(dist_L>1.5);
compare_probs_observer( MM, outliers, C )
    
%exportfig( gcf, [ 'outliers_boxplot.eps' ], 'Color', 'rgb' );

[jod, stats] = pw_scale_bootstrp( MM, bootstrap_samples );
%[jod, stats] = pw_scale_bootstrp( MM, bootstrap_samples, { 'use_parallel', 'never' } );

jod_offset = 1 - min(jod);
jod = jod + jod_offset;
stats.jod_high = stats.jod_high + jod_offset;
stats.jod_low = stats.jod_low + jod_offset;

error_bars = zeros(N, 2);
error_bars(:,1) = jod'-stats.jod_low';
error_bars(:,2) = stats.jod_high'-jod';

%html_change_figure_print_size( gcf, 22, 8 );
figure;
hold on;
barwitherr( error_bars, jod );
set( gca, 'XTickLabel', C );
grid on
xlabel( 'Conditions compared' );
ylabel( 'JOD' );
colormap(autumn)
hold off;

%exportfig( gcf, [ 'barplot.eps' ], 'Color', 'rgb' );
figure;
hold on;
pw_plot_ranking_triangles( jod, stats, C )
hold off;

%exportfig( gcf, [ 'triangles.eps' ], 'Color', 'rgb' );
