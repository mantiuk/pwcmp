% Another example showing how to execute scaling method with bootstrapping

% Add the scaling functions to the path if they are missing
if( ~exist( 'pw_scale', 'file' ) )
    addpath( '../' );
end

N_obs = 20; % Number of observers
N = 3; % Number of compared conditions

MM = zeros( N_obs, N^2 ); % We will store marices for all observers there

D_mean = [ 0  1  0;
           3  0  1;
           0  3  0 ];

% Let's generate comparison matrix for 10 observers
for kk=1:10
    
    M = max(D_mean + round( rand(size(D_mean))*3 ), 0 ); % Generate observer data with some noise    
    MM(kk,:) = M(:); % comparison matrix per observer "kk" is stored in the row "kk"
end
  
  
[Q, stats] = pw_scale_bootstrp( MM, 500 );

clf

eb_high = stats.jod_high - Q;
eb_low = Q - stats.jod_low;

errorbar( 1:N, Q, eb_low, eb_high );

xlim( [0.9 3.1] );


