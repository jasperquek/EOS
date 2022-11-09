function [] = plot_part_delta_vs_MLMCMC_obs_explicitRK3(full_delta, length_delta, starting_index_of_delta, trigger_condition, time_computation_mean, real_time_step, offline_data, E_MLMCMC_w_cov, ube_pop_var_time_to_erupt_allrepeats, ube_pop_var_final_tilt_all_repeats, L, rerun, part)

global length_effective_delta 

MCMC_time_level = 2;    % MCMC_time_step will be 0.5sec
factor_of_real_time = 4;
MCMC_time_step = factor_of_real_time*real_time_step*(2^(-MCMC_time_level));

len_para = length(E_MLMCMC_w_cov)-4;

[ observation_MLMCMC, time_of_erupt ] = MLMCMC_fun_observation_explicit_RK( offline_data, real_time_step, MCMC_time_level, E_MLMCMC_w_cov(1:len_para) );

% variance of time to eruption
% var_time_to_erupt = E_MLMCMC_w_cov(7)-(E_MLMCMC_w_cov(6)^2);
var_time_to_erupt = ube_pop_var_time_to_erupt_allrepeats;

% variance of max tilt
% var_max_tilt = E_MLMCMC_w_cov(9)-(E_MLMCMC_w_cov(8)^2);
var_max_tilt = ube_pop_var_final_tilt_all_repeats;

t_real = 0.5*( (1:length(full_delta)) - (starting_index_of_delta+length_delta-1) );
t_real_highlight = 0.5*( (-length_delta+1):0 );

delta = full_delta(starting_index_of_delta:(starting_index_of_delta+length_delta-1));

% Compute the difference of the mean of last 40% elements with mean of
% first 60% elements of delta
    [ diff_mean40_mean60 ] = func_diff_l40_f60( delta((end-length_effective_delta+1):end) );
    
% ---- Compute the approx time_index of full_observation that is equal to
% delta(end) ----
%     [approx_index_of_roots] = root_finder(delta, observation_MLMCMC, diff_mean40_mean60);
%     approx_index_of_roots = approx_index_of_roots + length_delta;
%     observation_MLMCMC = [zeros(length_delta,1);observation_MLMCMC];

    constants = {MCMC_time_step, real_time_step, length_delta, offline_data.constants{10}};
    [ Phi, time_to_erupt, observation_MLMCMC, approx_index_of_roots ] = MLMCMC_fun_potential2( constants, delta, observation_MLMCMC, time_of_erupt );
    
% hold on
t_MLMCMC = 0.5*( (1:length(observation_MLMCMC)) - approx_index_of_roots );
if approx_index_of_roots<=length_delta
    if diff_mean40_mean60 < trigger_condition
        t_MLMCMC_highlight = 0.5*( (1-approx_index_of_roots):0 );
        plot(t_real,full_delta,'k-',t_MLMCMC,observation_MLMCMC,'b-')
        hold on
        plot(t_real_highlight,delta,'r-',t_MLMCMC_highlight,observation_MLMCMC(1:approx_index_of_roots),'g-', 'LineWidth',1.5)
        title_var_time_erupt = ['Variance of time to eruption is ', num2str(var_time_to_erupt)];
%         title_var_max_tilt = ['Variance of max tilt is ', num2str(var_max_tilt)];
%         titlename = {title_var_time_erupt, title_var_max_tilt};
        titlename = title_var_time_erupt;
        title(titlename)
        xlabel('Time (sec)')
        ylabel('Tilt (nrad)')
    else
        % ---- Compute time to eruption from "start" ----
        time_to_erupt_ideal = (length(observation_MLMCMC) - approx_index_of_roots)*MCMC_time_step;
        time_left_to_run = time_to_erupt_ideal - time_computation_mean;
        
        t_MLMCMC_highlight = 0.5*( (1-approx_index_of_roots):0 );
        plot(t_real,full_delta,'k-',t_MLMCMC,observation_MLMCMC,'b-')
        hold on
        plot(t_real_highlight,delta,'r-',t_MLMCMC_highlight,observation_MLMCMC(1:approx_index_of_roots),'g-', 'LineWidth',1.5)
        title_var_time_erupt = ['Variance of time to eruption is ', num2str(var_time_to_erupt)];
%         title_var_max_tilt = ['Variance of max tilt is ', num2str(var_max_tilt)];
%         title_evacu = ['Run! Semeru will erupt in ', num2str(time_left_to_run), ' seconds!'];
        title_evacu = 'Sound off hypothetical evacuation alarm';
        title_forecast = ['Forecast eruption in ', num2str(time_left_to_run), ' seconds'];
%         titlename = {title_var_time_erupt, title_var_max_tilt, title_evacu};
        titlename = {title_var_time_erupt, title_evacu, title_forecast};
        title(titlename)
        xlabel('Time (sec)')
        ylabel('Tilt (nrad)')
    end
else
    if diff_mean40_mean60 < trigger_condition
        t_MLMCMC_highlight = 0.5*( (-length_delta+1):0 );
        plot(t_real,full_delta,'k-',t_MLMCMC,observation_MLMCMC,'b-')
        hold on
        plot(t_real_highlight,delta,'r-',t_MLMCMC_highlight,observation_MLMCMC(approx_index_of_roots-length_delta+1:approx_index_of_roots),'g-', 'LineWidth',1.5)
        title_var_time_erupt = ['Variance of time to eruption is ', num2str(var_time_to_erupt)];
%         title_var_max_tilt = ['Variance of max tilt is ', num2str(var_max_tilt)];
%         titlename = {title_var_time_erupt, title_var_max_tilt};
        titlename = title_var_time_erupt;
        title(titlename)
        xlabel('Time (sec)')
        ylabel('Tilt (nrad)')
    else
        % ---- Compute time to eruption from "start" ----
        time_to_erupt_ideal = (length(observation_MLMCMC) - approx_index_of_roots)*MCMC_time_step;
        time_left_to_run = time_to_erupt_ideal - time_computation_mean;
        
        t_MLMCMC_highlight = 0.5*( (-length_delta+1):0 );
        plot(t_real,full_delta,'k-',t_MLMCMC,observation_MLMCMC,'b-')
        hold on
        plot(t_real_highlight,delta,'r-',t_MLMCMC_highlight,observation_MLMCMC(approx_index_of_roots-length_delta+1:approx_index_of_roots),'g-', 'LineWidth',1.5)
        title_var_time_erupt = ['Variance of time to eruption is ', num2str(var_time_to_erupt)];
%         title_var_max_tilt = ['Variance of max tilt is ', num2str(var_max_tilt)];
%         title_evacu = ['Run! Semeru will erupt in ', num2str(time_left_to_run), ' seconds!'];
        title_evacu = 'Sound off hypothetical evacuation alarm';
        title_forecast = ['Forecast eruption in ', num2str(time_left_to_run), ' seconds'];
        titlename = {title_var_time_erupt, title_evacu, title_forecast};
        title(titlename)
        xlabel('Time (sec)')
        ylabel('Tilt (nrad)')
    end
end
hold off



filename = ['live_explicitRK_L', num2str(L),'_rerun', num2str(rerun),'_part',num2str(part)];

set(gcf, 'PaperPosition', [0 0 5.093 3.166]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5.093 3.166]); %Set the paper to have width 5 and height 5.
saveas(gcf, filename, 'pdf') %Save figure
saveas(gcf, filename) %Save figure

end