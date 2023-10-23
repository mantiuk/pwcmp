function [phi_rec, delta_rec, v_rec, OBSs, CONDs] = mos_scale( observer_id, condition_id, quality, options )
arguments
    observer_id (:,1)
    condition_id (:,1)
    quality (:,1) { mustBeNumeric }
    options.no_delta (1,1) = false
end
% The function implememts scaling of direct ranking data from the paper:
% Zhi Li, et al. "A Simple Model for Subject Behavior in Subjective
% Experiments" (https://arxiv.org/pdf/2004.02067)
%


if ~all(size(observer_id)==size(condition_id)) || ~all(size(condition_id)==size(quality))
    error( 'observer_id, condition_id and quality must be column vectors of the same size' );
end

[CONDs,~,cond_ind] = unique(condition_id);
N = numel(CONDs);
[OBSs,~,obs_ind] = unique(observer_id);
K = numel(OBSs);

phi_0 = ones(N,1) * mean(quality);
delta_0 = zeros(K,1);
log_v_0 = ones(K,1)*-2;

if options.no_delta

    par_0 = [phi_0' log_v_0'];

    par_rec = fminunc( @log_likelyhood_no_delta, par_0 );

    phi_rec = par_rec(1:N);
    v_rec = exp(par_rec((N+1):(N+K)));

else

    par_0 = [phi_0' delta_0' log_v_0'];
    %sigma_delta = (max(quality)-min(quality))/5;
    %alpha = 1; %1/(2*sigma_delta^2);
    par_rec = fminunc( @log_likelyhood, par_0 );

    phi_rec = par_rec(1:N);
    delta_rec = par_rec((N+1):(N+K));
    v_rec = exp(par_rec((N+K+1):(N+K*2)));

end

    function L = log_likelyhood( par )

        phi = par(1:N)';
        delta = par((N+1):(N+K))';
        log_v = par((N+K+1):(N+K*2))';

        L = sum( (quality-phi(cond_ind)-delta(obs_ind)).^2./(2*exp(log_v(obs_ind)).^2) + log_v(obs_ind) ); % + alpha*sum( delta.^2 );
    end

    function L = log_likelyhood_no_delta( par )

        phi = par(1:N)';
        log_v = par((N+1):(N+K))';

        L = sum( (quality-phi(cond_ind)).^2./(2*exp(log_v(obs_ind)).^2) + log_v(obs_ind) );
    end

end

