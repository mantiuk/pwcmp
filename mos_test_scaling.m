% Test scaling of MOS values

N = 70; % Number of assessed conditions
K = 15; % Number of observers
R = 5; % The number of repetitions

s = RandStream('mt19937ar', 'Seed', 10 );

phi = rand(s, N,1)*4+1;  % ground truth quality

sigma_delta = 0.5;
U_delta = randn(s, 1,K) * sigma_delta;
U_v = 1+rand(s, 1,K)*2;

U = phi + U_delta + U_v .* randn(s, N,K,R);

observer_id = repmat( 1:K, [N 1 R] );
condition_id = repmat( (1:N)', [1 K R] );

[phi_rec, delta_rec, v_rec, OBSs, CONDs] = mos_scale( observer_id(:), condition_id(:), U(:) );

%%
figure(1);
clf;
subplot(1,3,1);
rng = [min(phi) max(phi)];
plot( rng, rng, '--k' );
hold on;
scatter( phi, phi_rec );
Q_RMSE = sqrt(mean((phi(:)-phi_rec(:)).^2));
xlabel( 'Ground truth quality' )
ylabel( 'Scaled quality' )
title( sprintf( 'RMSE=%g', Q_RMSE) );

subplot(1,3,2);
rng = [min(U_delta) max(U_delta)];
plot( rng, rng, '--k' );
hold on;
scatter(U_delta, delta_rec);
xlabel( 'GT observer bias' );
ylabel( 'Scaled observer bias' );

subplot(1,3,3);
rng = [min(U_v) max(U_v)];
plot( rng, rng, '--k' );
hold on;
scatter(U_v, v_rec);
xlabel( 'GT observer error' );
xlabel( 'scaled observer error' );


ylabel( 'Observer''s std' );

