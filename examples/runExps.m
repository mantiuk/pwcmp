
res_file = 'simul_d_obs_withPrior.csv';

fh = fopen( res_file, 'w' );
fprintf( fh, 'N_obs, do_ties, full_design, bias, ci, d, mse\n' );


for dd=1:2,
    full_design = (dd==1);
    for tt=1 %(1:2)
        
        do_ties = (tt == 2);
        for oo=6:50
            N_obs = oo;
            
            q = [0 1 2 3 4];
            
            display( sprintf( 'N_obs = %d; ties = %d', N_obs, do_ties ) );
            [bias, ci, d, q_m, meas_mse] = simulate_exp_bootstrp( q, N_obs, do_ties, full_design );
            
            fprintf( fh, '%d, %d, %d, %g, %g, %g, %g\n', N_obs, do_ties, full_design, bias, ci, d, meas_mse );
            
        end
    end
end

fclose( fh );
