% Simple example showing how to execute scaling method

% Add the scaling functions to the path if they are missing
if( ~exist( 'pw_scale', 'file' ) )
    addpath( '../' );
end

% Comparison matrix
D = [ 0   25  25;
      75  0  25;
      75  75 0  ];
   
Q = pw_scale( D );

display( Q )
