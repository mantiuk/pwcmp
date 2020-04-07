function compare_probs_observer( MM, n_obs, C )
% Plots the probability of selecting any condition as better for all of the
% data and a selected observer. Ideally used with outlier analysis
%
% MM - responses for each observer in seperate rows.
% n_obs - number of observer to analyse
% C - cell array with the labels used for each plotted conditions
%

% Compute probabilities per observer

N = sqrt(size(MM,2));

for oo=1:size(MM,1)
    M = reshape(MM(oo,:),N,N);
    
    prob = M./(M + M');
    prob(isnan(prob)) = 0;
        
    s_columns = sum(prob,1);
    s_rows = sum(prob,2);
        
    data_boxplot_rows(oo,:) = s_rows'./(s_columns + s_rows');    
end

% Save the data from the outliers
o_rows = data_boxplot_rows(n_obs,:);    

% Remove the data from the outlier from the distribution
data_boxplot_rows(n_obs,:) = [];
%html_change_figure_print_size( gcf, 18, 8 );

figure;
hold on;
% Plot the data for that scene
boxplot(data_boxplot_rows,C)
plot([1:N],o_rows,'ok','Linewidth',2)

xlabel('Conditions compared')
ylabel('Probability of higher quality')
legend('Potential outlier')
%html_change_figure_print_size( gcf, 20, 8 );
grid on

hold off


end