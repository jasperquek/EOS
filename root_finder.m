function [approx_index_of_roots] = root_finder(delta, full_observation, diff_mean40_mean60)
    global trigger_condition
    global standard_dev
    if diff_mean40_mean60 < trigger_condition
%         full_observation_minus_delta_end = full_observation - mean( delta( (end-floor(0.1*length(delta))) :end) );
%         disp('no trigger')
        full_observation_minus_delta_end = full_observation - mean( delta( (end-10) :end) );
    else
%         disp('trigger')
        full_observation_minus_delta_end = full_observation - delta(end);
    end
    
    len_full_obs = length(full_observation);
    [ min_value, min_point ] = min( full_observation );
    no_of_roots = (full_observation_minus_delta_end(1)*full_observation_minus_delta_end(min_point)<0)...
                    +(full_observation_minus_delta_end(min_point)*full_observation_minus_delta_end(end)<0);
    
%     disp(no_of_roots)
    
    if no_of_roots == 1
        [approx_index_of_roots] = bisection_finder(full_observation_minus_delta_end, 1, len_full_obs );
    elseif no_of_roots == 2
%         if diff_mean40_mean60 < 0
        if delta(end)-min(delta) < 2*standard_dev
            [approx_index_of_roots] = bisection_finder(full_observation_minus_delta_end, 1, min_point );
        else
            [approx_index_of_roots] = bisection_finder(full_observation_minus_delta_end, min_point, len_full_obs );
        end
    elseif no_of_roots == 0
        [min_diff, approx_index_of_roots] = min(abs(full_observation_minus_delta_end));
    else
        error('There are more than 2 roots!')
    end
end