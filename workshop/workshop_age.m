%% Setting up the experiment
%clear all

% Set the path
set_path;

% Read aggragated results 
csv_file = readtable('aggregated_age.csv');

% Get observer ids 
observers = unique(csv_file.observer);

% Find all compared conditions
C = unique( cat( 1, csv_file.condition_1, csv_file.condition_2 ) );

% Find out the number of conditions in the experiment
N = numel(C);

% MM - KxC matrix with positive integers. Each row of the matrix contains
%      a comparison matrix for a single observer in 'flattenned' vectoried
%      format. That is, if a comparison matrix for observer k is M, MM 
%      should be initialized: 
%      MM(k,:) = M(:);
%      The number of columns C is equal to N*N where N is the number of
%      compared conditions. The number of rows K is equal the number of
%      observers.
MM = zeros(size(observers,1),N*N);

% For each observer
for ii = 1:size(observers,1)
    
    % Select only the rows correponding to a given observer
    ids_images = find(strcmp(csv_file.observer, observers{ii}));
    
    % Record results of each obsever to a matrix M
    M = zeros(N,N);
    
    % Go over all pairs of comparisons and record them into matrix M
    for jj = 1:numel(ids_images)
        kk = ids_images(jj);
        c1 = find( strcmp( csv_file.condition_1{kk}, C ) );
        c2 = find( strcmp( csv_file.condition_2{kk}, C ) );
        M(c1,c2) = M(c1,c2) + csv_file.chose_1(kk);
        M(c2,c1) = M(c2,c1) + csv_file.chose_2(kk);
    end
    
    % Flatten matrix M and put it into MM (as per description)
    MM(ii,:)=M(:)';   
end


%% Instructions

% You will need to complete parts of the code of this
% file. If you need any assistance, please raise your hand. 

%% Outlier detection and removal
% On this stage we perform a number of steps:      
%      1) Detect outliers in the matrix MM
%      2) Remove detected outliers from the matrix MM.
%      3) convert MM matrix into matrix D, (similar to matrix M, but  
%         containing no outliers).

%--------------------------------------------------------------------------

% Run the outlier analysis on matrix MM. Refer to the documentation of 
% pw_outlier_analysis function. Carefully study the input and the output of
% the function.

%Your code goes here ....
error( 'Code missing' );

%--------------------------------------------------------------------------

% Remove outliers from matrix MM.

% Note: you will need to study the output of pw_outlier_analysis to decide 
% on who are outliers. You can check ids of participants to find out who is 
% an outlier by selecting positions in the array of observer_ids.  
% e.g if outliers are in positions 1,5,7 type in oberver_ids{[1,5,7]}

% Hint: to remove rows 1,3,4 from the matrix K you need to type
% K([1,3,4],:) = [];

%Your code goes here ....
error( 'Code missing' );

%% Scaling with confidence intervals

% Once you have "clean" data without outliers, you can scale pairwise
% comparisons to obtain JND value for each tested condition (the photograph of
% a scientist). Please check the documentation of pw_scale_bootstrp
% function to figure out how to call it. The function can also estimate 
% confidence intervals using bootstrapping.

% Your code goes here ....
error( 'Code missing' );


%% Create a table with the results

age_scaled = table( C, jod, stats.jod_low, stats.jod_high, 'VariableNames', { 'name', 'jnd', 'jnd_low', 'jnd_high' } );

% Read true age data
true_age = readtable('true_ages.csv','Delimiter',',');

% Join the tables to add the true age
age_scaled = join( age_scaled, true_age );

% Save the scaled results to the CVS file so that it can be analysed with other
% software
writetable( age_scaled, 'age_scaled.csv', 'Delimiter', ',' );

%% Analysing the results

% Your final step is to plot the scaled results to obtain the plots similar 
% to the ones shown in the exercise description. You can use either Matlab
% (recommended), or you favourite plotting software. 

% Please prepare:
% 1) a scatter plot of "true age" vs. JND of age difference;
% 2) add confidence intervals to the plot;
% 3) compute the correlation between the true age and JND values.
%
% The Matlab functions that could be useful for this task: scatter, errorbar,
% text, corr
% 
% Given the plot, find for which photograph the age was most often
% confused. 

%Your code goes here ....





