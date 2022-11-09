% @@@@@@@@@@ --- Function to compute potential --- @@@@@@@@@@@@@@@@@@@@@@@@
function [ Phi, time_to_erupt, full_observation_pad_zeros, better_root ] = MLMCMC_fun_potential2( constants, delta, full_observation, time_of_erupt )
%UNTITLED Summary of this function goes here
%   constants - {MCMC_time_step, real_time_step, no_of_obs, obs_cov_mat }
%   delta - part of the full time series of tilt, since in practice, we
%           only have a partial information of the full explosion

global length_effective_delta

MCMC_time_step = constants{1};
real_time_step = constants{2};
no_of_obs = constants{3};
obs_cov_mat = constants{4};

% Compute the difference of the mean of last 40% elements with mean of
% first 60% elements of delta, of the last 42sec of delta
    [ diff_mean40_mean60 ] = func_diff_l40_f60( delta((end-length_effective_delta+1):end) );
    
    
% ---- Pad full observations with zeros -----------
    full_observation_pad_zeros = [ zeros( ceil( (no_of_obs-1)*real_time_step/MCMC_time_step ), 1); full_observation ];


% ---- Compute the approx time_index of full_observation_pad_zeros that is equal to delta(end) ----
    [approx_index_of_roots] = root_finder(delta, full_observation_pad_zeros, diff_mean40_mean60);
    if approx_index_of_roots == 1
        % points to the first element of full_observation, before padding 
        approx_index_of_roots = ceil( (no_of_obs-1)*real_time_step/MCMC_time_step ) + 1;
    end


% ------ (approx_index_of_roots-1)*MCMC_time_step is always MORE THAN OR EQUAL  to
% (no_of_obs-1)*real_time_step --------------------
%
    % % ---- Check if time before approx_index_of_roots is less than duration of
    % % delta ----
    %     if (approx_index_of_roots-1)*MCMC_time_step <= (no_of_obs-1)*real_time_step
    %         observation = full_observation(1:approx_index_of_roots);
    %     else
    %         observation = full_observation( (approx_index_of_roots-floor((no_of_obs-1)*real_time_step/MCMC_time_step)) : approx_index_of_roots );
    %     end
    
    % If MCMC_time_step <= real_time_step, then
    %   (no_of_obs-1)*real_time_step/MCMC_time_step is an integer, that is
    %       length of time of observation will be equal to length of time of delta
    % else MCMC_time_step > real_time_step
    %   (no_of_obs-1)*real_time_step/MCMC_time_step might have a remainder which
    %   we do not want. 
    %       i.e. length of time of observation will be less
    %       than length of time of delta
    observation = full_observation_pad_zeros( (approx_index_of_roots-floor((no_of_obs-1)*real_time_step/MCMC_time_step)) : approx_index_of_roots );
    no_syn_obs = length(observation);

    
if MCMC_time_step>real_time_step    %there's more real life data than synthetic data
    if (no_of_obs-1)*real_time_step < (no_syn_obs-1)*MCMC_time_step     %------------->this case will not happen
        disp('I am wrong, real time is less than mcmc time')
        %delta_at_MCMC_time_step = flip( delta(no_of_obs:(-(MCMC_time_step/real_time_step)):1) );    %flip to maintain chronological order
        %delta_padded = [ zeros( (no_syn_obs-length(delta_at_MCMC_time_step)), 1); delta_at_MCMC_time_step];     %in chronological order, same legnth as observation
        %Phi = 0.5*( (delta_padded-observation)' * (delta_padded-observation) )/obs_cov_mat;
        %-----
        
    elseif (no_of_obs-1)*real_time_step > (no_syn_obs-1)*MCMC_time_step 	%------------->this might be the case
        [Phi, better_root] = neighbor_search_better_Phi_MgR(approx_index_of_roots, MCMC_time_step, full_observation_pad_zeros, real_time_step, no_of_obs, delta, obs_cov_mat);
        %-----
    else %------------->this might be the case
        [Phi, better_root] = neighbor_search_better_Phi_MgR(approx_index_of_roots, MCMC_time_step, full_observation_pad_zeros, real_time_step, no_of_obs, delta, obs_cov_mat);
        %-----
    end
elseif MCMC_time_step<real_time_step    %there's more synthetic life data than real data
    if (no_of_obs-1)*real_time_step < (no_syn_obs-1)*MCMC_time_step     %------------->this will not happen
        disp('I am wrong, real time is less than mcmc time')
        %observation_at_real_time_step = flip( observation(  no_syn_obs : -(real_time_step/MCMC_time_step) : 1  ) );  %flip to maintain chronological order
        %delta_padded = [ zeros( (length(observation_at_real_time_step) - no_of_obs), 1); delta ];
        %Phi = 0.5*( (delta_padded-observation_at_real_time_step)' * (delta_padded-observation_at_real_time_step) )/obs_cov_mat;
        %-----
        
    elseif (no_of_obs-1)*real_time_step > (no_syn_obs-1)*MCMC_time_step     %------------->this will not happen
        disp('I am wrong, real time is more than mcmc time')
        %observation_at_real_time_step = flip( observation(  no_syn_obs : -(real_time_step/MCMC_time_step) : 1  ) );  %flip to maintain chronological order
        %observation_padded = [ zeros( (no_of_obs-length(observation_at_real_time_step)), 1 ); observation_at_real_time_step];
        %Phi = 0.5*( (delta-observation_padded)' * (delta-observation_padded) )/obs_cov_mat;
        %-----
        
    else %------------->this is always the case
        [Phi, better_root] = neighbor_search_better_Phi_MleqR(approx_index_of_roots, MCMC_time_step, full_observation_pad_zeros, real_time_step, no_of_obs, delta, obs_cov_mat);
        %-----
    end
else %MCMC_time_step == real_time_step
    if (no_of_obs < no_syn_obs)     %------------->this will not happen
        disp('I am wrong, real time is less than mcmc time')
%         Phi = 0.5*((delta - observation( (no_syn_obs-no_of_obs+1):no_syn_obs) )')*(obs_cov_mat\(delta - observation( (no_syn_obs-no_of_obs+1):no_syn_obs) ) );
        %delta_padded = [ zeros( (no_syn_obs-no_of_obs), 1 ); delta];
        %Phi = 0.5*( (delta_padded-observation)' * (delta_padded-observation) )/obs_cov_mat;
    elseif (no_of_obs > no_syn_obs)     %------------->this will not happen
        disp('I am wrong, real time is less than mcmc time')
%         Phi = 0.5*((delta( (no_of_obs-no_syn_obs+1):no_of_obs ) - observation)')*(obs_cov_mat( (no_of_obs-no_syn_obs+1):no_of_obs, (no_of_obs-no_syn_obs+1):no_of_obs )\(delta( (no_of_obs-no_syn_obs+1):no_of_obs ) - observation ) );
        %observation_padded = [ zeros( (no_of_obs-no_syn_obs), 1 ); observation];
        %Phi = 0.5*( (delta-observation_padded)' * (delta-observation_padded) )/obs_cov_mat;
    else    %------------->this is always the case
        [Phi, better_root] = neighbor_search_better_Phi_MleqR(approx_index_of_roots, MCMC_time_step, full_observation_pad_zeros, real_time_step, no_of_obs, delta, obs_cov_mat);
    end
end  

% ----- Compute the time to eruption -------
    time_to_erupt = time_of_erupt - (better_root-1 - ceil( (no_of_obs-1)*real_time_step/MCMC_time_step))*MCMC_time_step;

end
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [better_Phi, better_root] = neighbor_search_better_Phi_MgR(approx_index_of_roots, MCMC_time_step, full_observation_pad_zeros, real_time_step, no_of_obs, delta, obs_cov_mat)
% When MCMC_time_step is greater than real_time_step

    % Radius of time to search
    radius_of_time = 100;
    
    % default index of full_observation_pad_zeros to start with
    index_wrt_approx_root = (approx_index_of_roots-floor((no_of_obs-1)*real_time_step/MCMC_time_step)) : approx_index_of_roots;
    
    % Length of time available to search on the left
    avail_time_on_left = ( index_wrt_approx_root(1) - 1 )*MCMC_time_step;
    avail_time_on_left = (avail_time_on_left <= radius_of_time)*avail_time_on_left + (avail_time_on_left > radius_of_time)*radius_of_time;
    
    % Length of time available to search on the right
    avail_time_on_right = ( length(full_observation_pad_zeros) - approx_index_of_roots )*MCMC_time_step;
    avail_time_on_right = (avail_time_on_right <= radius_of_time)*avail_time_on_right + (avail_time_on_right > radius_of_time)*radius_of_time;
    
    % Number of index avaiable to move to the left and right
    avail_index_left = avail_time_on_left/MCMC_time_step;
    avail_index_right = avail_time_on_right/MCMC_time_step;
    
    % delta at MCMC time steps
    delta_at_MCMC_time_step = flip( delta(no_of_obs:(-(MCMC_time_step/real_time_step)):1) );
    
    % Collection of Candidates tilts (observation), centered at
    % approx_in_of_roots
    candidate_potentials = zeros( length(delta_at_MCMC_time_step), (avail_index_left+avail_index_right+1) );
    for i = 1:(avail_index_left+avail_index_right+1)
        candidate_potentials(:,i) = delta_at_MCMC_time_step - full_observation_pad_zeros( index_wrt_approx_root + (-avail_index_left+(i-1)) );
    end
    candidate_potentials = sum(candidate_potentials.^2)/(2*obs_cov_mat);
    
    % Choose the minimum potential 
    [better_Phi, better_root] = min(candidate_potentials);
    better_root = approx_index_of_roots + (-avail_index_left+(better_root-1));
end

function [better_Phi, better_root] = neighbor_search_better_Phi_MleqR(approx_index_of_roots, MCMC_time_step, full_observation_pad_zeros, real_time_step, no_of_obs, delta, obs_cov_mat)
% When MCMC_time_step is greater than real_time_step

    % Radius of time to search
    radius_of_time = 100;
    
    % default index of full_observation_pad_zeros to start with
    index_wrt_approx_root = (approx_index_of_roots-floor((no_of_obs-1)*real_time_step/MCMC_time_step)): (real_time_step/MCMC_time_step): approx_index_of_roots;
    
    % Length of time available to search on the left
    avail_time_on_left = ( index_wrt_approx_root(1) - 1 )*MCMC_time_step;
    avail_time_on_left = (avail_time_on_left <= radius_of_time)*avail_time_on_left + (avail_time_on_left > radius_of_time)*radius_of_time;
    
    % Length of time available to search on the right
    avail_time_on_right = ( length(full_observation_pad_zeros) - approx_index_of_roots )*MCMC_time_step;
    avail_time_on_right = (avail_time_on_right <= radius_of_time)*avail_time_on_right + (avail_time_on_right > radius_of_time)*radius_of_time;
    
    % Number of index avaiable to move to the left and right, with jumping
    % by real_time_step
    avail_jump_index_left = avail_time_on_left/real_time_step;
    avail_jump_index_right = avail_time_on_right/real_time_step;
    
    % Collection of Candidates tilts (observation), centered at
    % approx_in_of_roots
    candidate_potentials = zeros( length(delta), (avail_jump_index_left+avail_jump_index_right+1) );
    for i = 1:(avail_jump_index_left+avail_jump_index_right+1)
        candidate_potentials(:,i) = delta - full_observation_pad_zeros( index_wrt_approx_root + ( (-avail_jump_index_left+(i-1))*(real_time_step/MCMC_time_step) ) );
    end
    candidate_potentials = sum(candidate_potentials.^2)/(2*obs_cov_mat);
    
    % Choose the minimum potential 
    [better_Phi , better_root] = min(candidate_potentials);
    better_root = approx_index_of_roots + (-avail_jump_index_left+(better_root-1))*(real_time_step/MCMC_time_step);
end