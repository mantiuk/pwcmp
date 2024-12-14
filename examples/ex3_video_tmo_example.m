% Analyse the data in the video_TMO example: outlier detection, scaling and
% statistical analysis

if( ~exist( 'pw_scale', 'file' ) )
    addpath( genpath('../') );
end

bootstrap_samples = 500;

D = readtable( 'ex3_tmo_cmp_data.csv' );

C = unique( cat( 1, D.condition_A, D.condition_B ) ); % all conditions


figure(1);
[L,dist_L, OBSs] = pw_outlier_analysis_table( D, { 'scene' }, { 'condition_A', 'condition_B' }, 'observer', 'is_A_selected' );

% Remove observers for which dist_L>1.5
fprintf( 1, 'Detected outliers:\n' );
fprintf( 1, '\t%s\n', OBSs{dist_L>1.5} );
D_no_outlier = D(~ismember(D.observer, OBSs(dist_L>1.5)),:);

%% Scale the pairwise comparisons (excluding outliers)

[R, Rs] = pw_scale_table(D_no_outlier, 'scene', { 'condition_A', 'condition_B' }, ...
    'observer', 'is_A_selected', 'bootstrap_samples', bootstrap_samples, 'do_all', true );

SCENEs = unique(R.scene);
JODs = zeros(length(SCENEs), length(C));
JODs_eb = zeros(length(SCENEs), length(C), 2); % error bars
for ss=1:length(SCENEs)
    for cc=1:length(C)
        sel = strcmp(R.scene,SCENEs{ss}) & strcmp(R.condition,C{cc});
        JODs(ss,cc) = R.jod(sel);
        JODs_eb(ss,cc,1) = R.jod(sel) - R.jod_low(sel);
        JODs_eb(ss,cc,2) = R.jod_high(sel) - R.jod(sel);
    end
end

% Add offset to have only positive JODs
jod_offset = - min(JODs) + 1;
JODs = JODs + jod_offset;


figure(2);
clf;
hold on;
barwitherr(JODs_eb, JODs);
xticklabels( SCENEs );
grid on
ylabel( 'Quality [JOD]' );
colormap(autumn)
hold off;
legend( C, 'Interpreter', 'none', 'Location', 'eastoutside' );

%% Perform statistical test for `all` group
figure(3);
pw_plot_ranking_triangles( Rs{1}.jod, Rs{1}.stats, C )

