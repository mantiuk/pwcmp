function pw_plot_ranking_triangles( jnd, stats, CONDs )
% Visualize significant differences using triangle plots
% 
% pw_plot_ranking_triangles( jnd, stats, CONDs )
%
% jnd - JND values, produced by pw_scale_bootstrp
% stats - stats structure, produced by pw_scale_bootstrp. The only relevant
%         element of that structure is covariance matrix jnd_cov.
% CONDs - cell array with the labels used for each plotted conditions. Must be the 
%         same size as "jnd".


% Construct a matrix with statistically significant differences between
% pairs of conditions
[p, D] = pw_significance_matrix( jnd, stats );
ssd = p<0.025; % Two-tailed test

scs = jnd;
cond_count = length(jnd);

[~, rank] = sort( scs, 'ascend' );

xy = zeros(cond_count,2);

for k=1:cond_count
    j = rank(k);
    xy(j,1) = scs(j);
    xy(j,2) = 1-mod( k, 2 );

    plot( xy(j,1), xy(j,2), 'or' );
    hold on
    offset = -(-1)^(mod( k, 2 )+1) * 0.1;
    text( xy(j,1), xy(j,2)+offset, CONDs{j}, 'HorizontalAlignment', 'center' );

    if( k>1 )
        draw_edge( j, rank(k-1) );
    end
    if( k > 2 )
        draw_edge( j, rank(k-2) );
    end
    
end

hold off
ylim( [-0.5 1.5 ] );
set( gca, 'YTick', [] );

    function draw_edge( v1, v2 )
%        if( edges_done(v1,v2) ), return, end;
        
        str = sprintf( '%.2f', abs(D(v1,v2)) );
        power = 1;
        if( ssd(v1,v2) )
            color = '-b';
        else
            color = '--r';
            
            % power analysis
%             sigma = sqrt(aDS.(['var_' score_column])(v1) + aDS.(['var_' score_column])(v2));
%             delta = aDS.(['mean_' score_column])(v1) - aDS.(['mean_' score_column])(v2);
%             nel = aDS.GroupCount(v1) + aDS.GroupCount(v2);
%             if( sigma ~= 0 )
%                 power = sampsizepwr( 't', [0 sigma], abs(delta), [], nel );
%             end
%             
%             if( power < 0.8 && sigma ~= 0 && delta>0 )                
%                 n_d = sampsizepwr( 't', [0 sigma], abs(delta), 0.8 );
%                 str = [str sprintf( ' [P(%d)=%.2g n_d=%d]', nel, power, n_d )];
%             else
%                 str = [str sprintf( ' [P(%d)=%.2g]', nel, power )];
%             end
        end
        
        plot( [xy(v1,1) xy(v2,1)], [xy(v1,2) xy(v2,2)], color );
        px = mean([xy(v1,1) xy(v2,1)]);
        py = mean([xy(v1,2) xy(v2,2)]);
                
        text( px, py, str, ...
            'HorizontalAlignment', 'center', 'Background', 'white' );
        
%        edges_done(v1,v2) = true;
    end


end


function str = all2str( v )

if( isnumeric( v ) )
    str = num2str( v );
else
    str = char( v );
end

end
