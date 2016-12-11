% Simple example showing how to execute scaling method

% Add the scaling functions to the path if they are missing
if( ~exist( 'pw_scale', 'file' ) )
    addpath( '../' );
end

% Comparison matrix. This is incomplete design in which 1-st conditon was
% compared with the 2nd and 2nd with the 3rd. In both cases 75 observers
% seleted one condition and 25 the other.
D = [ 0   25  0;
      75  0  25;
      0   75 0  ];
   
Q = pw_scale( D );

display( Q )
