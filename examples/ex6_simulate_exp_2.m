% Simulate experiment before it is run

m = 3;
m_size = 4;

do_neighboring_levels = false;

q = [0 1:4 (1:4)*1.1 (1:4)*1.5];

C_mat = false(length(q));
for kk=1:m
    ind = 2+(kk-1)*m_size;
    C_mat(1,ind) = true; % Comparison with the reference
    for ll=1:(m_size-1)
        ind = 1+(kk-1)*m_size+ll;
        C_mat(ind,ind+1) = true; % Within the method
    end
end


for ll=1:m_size
    for kk=1:m
        for ii=(kk+1):m
            ind1 = 1+(kk-1)*m_size+ll;
            ind2 = 1+(ii-1)*m_size+ll;
            C_mat(ind1,ind2) = true;   % Between the methods
            
            if( do_neighboring_levels )
                if( ll>1 )
                    C_mat(ind1,ind2-1) = true;   
                end
                if( ll<m_size )
                    C_mat(ind1,ind2+1) = true;   
                end
            end
        end
    end
end

      
q_s = pw_simulate_exp( 1000, q, 40, C_mat );


COLORs = lines(m);

q_mean = mean(q_s);
q_std = std(q_s);

for kk=1:m
    ind = (kk-1)*m_size+2;
    rng = ind:(ind+m_size-1);
    plot( 0:m_size, [0 q(rng)], '-o', 'Color', COLORs(kk,:) );
    hold on
    errorbar( 0:m_size, [0 q_mean(rng)], [0 q_std(rng)], '--x', 'Color', COLORs(kk,:) );
end
hold off
xlabel( 'Distortion level' );
ylabel( 'Quality' );

grid on



display( sprintf( 'Mean std: %g', mean(q_std(2:end).^2).^0.5 ) );

