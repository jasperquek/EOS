% ----------- Generating Synthetic Reference Observation-------------
    load('Offline_variables')
    load('Online_variables')
    
    % ------------- MLMCMC parameters ----------------
    real_time_step = 0.5;
    MCMC_time_level = 2;

    % ---- Initializing constants for properties of online parameters ----
    online_parameters = online_data;
    %------------------------------------------------

    [ tilt_total ] = MLMCMC_fun_observation_explicit_RK( offline_data, real_time_step, MCMC_time_level, online_parameters );

    % Noise of tilt at any future time is assume to distribute identically when
    % t=0
    standard_deviation = tilt_total(1);

    % Generate delta (noisy observation)
    full_delta = tilt_total + standard_deviation*randn(size(tilt_total));
    
    % pad synthetic reference observation with only noises from 1st 750sec
    full_delta = [ standard_deviation*randn(1501,1); full_delta];

    % Covariance matrix
    % obs_cov_mat = diag( (standard_deviation^2) * ones(1,length(tilt_total)) );
    obs_cov_mat = standard_deviation^2;

    filename = 'fulldelta_and_obscovmat.mat';
    save(filename,'full_delta','obs_cov_mat')