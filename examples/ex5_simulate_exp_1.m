% Simulate multiple experiment runs before of after the actual experiment
% is run
%
% This is can help to determine how many observers are needed and which
% pairs of conditions should be compared in the incomplete design. 

q = [0 1 2 3];

C_mat = [ 0 1 0 0;
          0 0 1 0;
          0 0 0 1;
          0 0 0 0 ];
      
q_s = pw_simulate_exp( 100, q, 20, C_mat );

clf

q_mean = mean(q_s);
q_std = std(q_s);

plot( [0 max(q)], [0 max(q)], '--k' );
hold on
errorbar( q, q_mean, q_std );
hold off
xlabel( 'True values' );
ylabel( 'Estimated values' );

grid on
