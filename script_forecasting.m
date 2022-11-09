%%
% ---- Trigger condition, in terms of standard deviation of noise, in noise free condition ----
    trigger_factor = 2;
    
% ----- Registering physical computing resource -----
    % if no_of_cores > 1, code will repeat the algorithm = no_of_cores and
    % take the average.
    no_of_cores = 1;

%%
% ----- loading offline_data and synthetic reference observation -------------
    load('fulldelta_and_obscovmat')
    load('Offline_variables')
    offline_data.constants{end+1} = obs_cov_mat;

    global standard_dev
    standard_dev = sqrt(obs_cov_mat);

% ----- Defining length of observation in practice and the time step for each update --------------
    length_full_delta = length(full_delta);
    length_delta = 1501;
    
    % length_effective_delta is the last 42sec of the delta that we will
    % use to check if it triggers the alarm
    global length_effective_delta 
    length_effective_delta = 84;  %each time step is 0.5sec, so this is 42sec duration

    starting_index_of_delta = 1;
    delta_step = 50;    %ideally, this should be set to be equal to computation time
    
% ----- Trigger factor, in terms of standard deviation, 
% ----- of the difference of the mean(first 60% of delta)-mean(last 40%)
% ----- in noisy condition. i.e. if diff_l40-f60 > trigger_factor_noisy 
% ----- then alarm will go off -----
    length_first_60_delta = ceil(0.6*length_effective_delta);
    global trigger_condition 
    trigger_condition = ( trigger_factor * sqrt(1/length_first_60_delta + 1/(length_effective_delta - length_first_60_delta)) )*sqrt(obs_cov_mat);


% ----- MLMCMC parameters -----
    real_time_step = 0.5;
    burn = 100;
    L = 4;
    rate = 1;
    starting_level = 1;
    alpha = 0;
    repeat = no_of_cores;


%---- File naming function -----
    rerun = 1;
    part = 1;
    filename = ['MLMCMC_explicitRK_a', num2str(alpha), '_L', num2str(L),'_rerun', num2str(rerun),'_part',num2str(part),'.mat'];
    while exist(filename,'file')==2
        rerun = rerun + 1;
        filename = ['MLMCMC_explicitRK_a', num2str(alpha), '_L', num2str(L),'_rerun', num2str(rerun),'_part',num2str(part),'.mat'];
    end
% ------------------------------



%% ----- Main body -------------
    online_parameters = gen_proposal(offline_data);
    if repeat == 1  % using only 1 computer core
        while (starting_index_of_delta+length_delta-1)<=length_full_delta
            % ------ Updating delta ----------------
                delta = full_delta(starting_index_of_delta:(starting_index_of_delta+length_delta-1));
                if length(delta) ~= length_delta
                    err_msg = ['Length of delta is ', num2str(length(delta)),'. It should be ',num2str(length_delta)];
                    error(err_msg)
                end

            %---- File naming function -----
                while exist(filename,'file')==2
                    part = part + 1;
                    filename = ['MLMCMC_explicitRK_a', num2str(alpha), '_L', num2str(L),'_rerun', num2str(rerun),'_part',num2str(part),'.mat'];
                end

            % ---- Main MLMCMC ---------
                E_MLMCMC = zeros(repeat, length(online_parameters)+4);
                time = zeros(repeat,1);
                ube_pop_var_time_to_erupt = zeros(repeat,1);
                ube_pop_var_final_tilt = zeros(repeat,1);
                for i = 1:repeat
                    tic
                    [E_MLMCMC(i,:), ube_pop_var_time_to_erupt(i), ube_pop_var_final_tilt(i)] = MLMCMC_Semeru_uniform_prior_explicitRK_1burn_col_timetilt( offline_data, online_parameters, delta, L, rate, starting_level, real_time_step, alpha, burn );
                    time(i) = toc;
                end

            % ----- saving variables -----
                E_MLMCMC_mean = E_MLMCMC;
                ube_pop_var_time_to_erupt_allrepeats = mean(ube_pop_var_time_to_erupt);
                ube_pop_var_final_tilt_all_repeats = mean(ube_pop_var_final_tilt);
                time_computation_mean = mean(time);
                save(filename,'E_MLMCMC_mean','E_MLMCMC','ube_pop_var_time_to_erupt_allrepeats','ube_pop_var_final_tilt_all_repeats','rate','L','alpha','time','delta','starting_index_of_delta','trigger_condition','time_computation_mean', 'real_time_step', 'offline_data', 'rerun', 'part')

            % ----- plotting -----
                plot_part_delta_vs_MLMCMC_obs_explicitRK3(full_delta, length_delta, starting_index_of_delta, trigger_condition, time_computation_mean, real_time_step, offline_data, E_MLMCMC_mean, ube_pop_var_time_to_erupt_allrepeats, ube_pop_var_final_tilt_all_repeats, L, rerun, part)

            % ---- Initializing constants for properties of online parameters for next set of delta ----
                online_parameters = E_MLMCMC_mean(1:length(online_parameters));
                starting_index_of_delta = starting_index_of_delta + delta_step;
        end

        if starting_index_of_delta-delta_step+length_delta-1 < length_full_delta
            starting_index_of_delta = length_full_delta - length_delta + 1;
            % ------ Updating delta ----------------
                delta = full_delta(starting_index_of_delta:(starting_index_of_delta+length_delta-1));
                if length(delta) ~= length_delta
                    err_msg = ['Length of delta is ', num2str(length(delta)),'. It should be ',num2str(length_delta)];
                    error(err_msg)
                end

            % ----- File naming function ----------
                while exist(filename,'file')==2
                    part = part + 1;
                    filename = ['MLMCMC_explicitRK_a', num2str(alpha), '_L', num2str(L),'_rerun', num2str(rerun),'_part',num2str(part),'.mat'];
                end

            % ------ Main MLMCMC ---------
                E_MLMCMC = zeros(repeat, length(online_parameters)+4);
                ube_pop_var_time_to_erupt = zeros(repeat,1);
                ube_pop_var_final_tilt = zeros(repeat,1);
                time = zeros(repeat,1);
                for i = 1:repeat
                    tic
                    [E_MLMCMC(i,:), ube_pop_var_time_to_erupt(i), ube_pop_var_final_tilt(i)] = MLMCMC_Semeru_uniform_prior_explicitRK_1burn_col_timetilt( offline_data, online_parameters, delta, L, rate, starting_level, real_time_step, alpha, burn );
                    time(i) = toc
                end

            % ----- saving variables -----
                E_MLMCMC_mean = E_MLMCMC;
                ube_pop_var_time_to_erupt_allrepeats = mean(ube_pop_var_time_to_erupt);
                ube_pop_var_final_tilt_all_repeats = mean(ube_pop_var_final_tilt);
                time_computation_mean = mean(time);
                save(filename,'E_MLMCMC_mean','E_MLMCMC','ube_pop_var_time_to_erupt_allrepeats','ube_pop_var_final_tilt_all_repeats','rate','L','alpha','time','delta','starting_index_of_delta','trigger_condition','time_computation_mean', 'real_time_step', 'offline_data', 'rerun', 'part')

            % ----- plotting -----
                plot_part_delta_vs_MLMCMC_obs_explicitRK3(full_delta, length_delta, starting_index_of_delta, trigger_condition, time_computation_mean, real_time_step, offline_data, E_MLMCMC_mean, ube_pop_var_time_to_erupt_allrepeats, ube_pop_var_final_tilt_all_repeats, L, rerun, part)
        end
    else  % using more than 1 computer core
        while (starting_index_of_delta+length_delta-1)<=length_full_delta
            % ------ Updating delta ----------------
                delta = full_delta(starting_index_of_delta:(starting_index_of_delta+length_delta-1));
                if length(delta) ~= length_delta
                    err_msg = ['Length of delta is ', num2str(length(delta)),'. It should be ',num2str(length_delta)];
                    error(err_msg)
                end

            %---- File naming function -----
                while exist(filename,'file')==2
                    part = part + 1;
                    filename = ['MLMCMC_explicitRK_a', num2str(alpha), '_L', num2str(L),'_rerun', num2str(rerun),'_part',num2str(part),'.mat'];
                end

            % ---- Main MLMCMC ---------
                E_MLMCMC = zeros(repeat, length(online_parameters)+4);
                time = zeros(repeat,1);
                ube_pop_var_time_to_erupt = zeros(repeat,1);
                ube_pop_var_final_tilt = zeros(repeat,1);
                parfor i = 1:repeat
                    tic
                    [E_MLMCMC(i,:), ube_pop_var_time_to_erupt(i), ube_pop_var_final_tilt(i)] = MLMCMC_Semeru_uniform_prior_explicitRK_1burn_col_timetilt( offline_data, online_parameters, delta, L, rate, starting_level, real_time_step, alpha, burn );
                    time(i) = toc
                end

            % ----- saving variables -----
                E_MLMCMC_mean = mean(E_MLMCMC);
                ube_pop_var_time_to_erupt_allrepeats = mean(ube_pop_var_time_to_erupt);
                ube_pop_var_final_tilt_all_repeats = mean(ube_pop_var_final_tilt);
                time_computation_mean = mean(time);
                save(filename,'E_MLMCMC_mean','E_MLMCMC','ube_pop_var_time_to_erupt_allrepeats','ube_pop_var_final_tilt_all_repeats','rate','L','alpha','time','delta','starting_index_of_delta','trigger_condition','time_computation_mean', 'real_time_step', 'offline_data', 'rerun', 'part')

            % ----- plotting -----
                plot_part_delta_vs_MLMCMC_obs_explicitRK3(full_delta, length_delta, starting_index_of_delta, trigger_condition, time_computation_mean, real_time_step, offline_data, E_MLMCMC_mean, ube_pop_var_time_to_erupt_allrepeats, ube_pop_var_final_tilt_all_repeats, L, rerun, part)

            % ---- Initializing constants for properties of online parameters for next set of delta ----
                online_parameters = E_MLMCMC_mean(1:length(online_parameters));
                starting_index_of_delta = starting_index_of_delta + delta_step;
        end

        if starting_index_of_delta-delta_step+length_delta-1 < length_full_delta
            starting_index_of_delta = length_full_delta - length_delta + 1;
            % ------ Updating delta ----------------
                delta = full_delta(starting_index_of_delta:(starting_index_of_delta+length_delta-1));
                if length(delta) ~= length_delta
                    err_msg = ['Length of delta is ', num2str(length(delta)),'. It should be ',num2str(length_delta)];
                    error(err_msg)
                end

            % ----- File naming function ----------
                while exist(filename,'file')==2
                    part = part + 1;
                    filename = ['MLMCMC_explicitRK_a', num2str(alpha), '_L', num2str(L),'_rerun', num2str(rerun),'_part',num2str(part),'.mat'];
                end

            % ------ Main MLMCMC ---------
                E_MLMCMC = zeros(repeat, length(online_parameters)+4);
                ube_pop_var_time_to_erupt = zeros(repeat,1);
                ube_pop_var_final_tilt = zeros(repeat,1);
                time = zeros(repeat,1);
                parfor i = 1:repeat
                    tic
                    [E_MLMCMC(i,:), ube_pop_var_time_to_erupt(i), ube_pop_var_final_tilt(i)] = MLMCMC_Semeru_uniform_prior_explicitRK_1burn_col_timetilt( offline_data, online_parameters, delta, L, rate, starting_level, real_time_step, alpha, burn );
                    time(i) = toc
                end

            % ----- saving variables -----
                E_MLMCMC_mean = mean(E_MLMCMC);
                ube_pop_var_time_to_erupt_allrepeats = mean(ube_pop_var_time_to_erupt);
                ube_pop_var_final_tilt_all_repeats = mean(ube_pop_var_final_tilt);
                time_computation_mean = mean(time);
                save(filename,'E_MLMCMC_mean','E_MLMCMC','ube_pop_var_time_to_erupt_allrepeats','ube_pop_var_final_tilt_all_repeats','rate','L','alpha','time','delta','starting_index_of_delta','trigger_condition','time_computation_mean', 'real_time_step', 'offline_data', 'rerun', 'part')

            % ----- plotting -----
                plot_part_delta_vs_MLMCMC_obs_explicitRK3(full_delta, length_delta, starting_index_of_delta, trigger_condition, time_computation_mean, real_time_step, offline_data, E_MLMCMC_mean, ube_pop_var_time_to_erupt_allrepeats, ube_pop_var_final_tilt_all_repeats, L, rerun, part)
        end
    end