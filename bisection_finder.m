function [root] = bisection_finder(full_observation_minus_delta_end, start_index, end_index )
    if (full_observation_minus_delta_end(start_index)*full_observation_minus_delta_end(end_index)>0)
        err_msg = 'Starting and ending index have the same sign';
        error(err_msg)
    end
    root = ceil( (start_index + end_index)/2 );
    
    counter = 1;
    iteration_limit = 2*log2(end_index-start_index);
    while not(end_index-start_index==1)
        % ---- Bisection method ---------
        if (full_observation_minus_delta_end(start_index)*full_observation_minus_delta_end(root)>0)
            start_index = root;
        else
            end_index = root;
        end
        root = ceil( (start_index + end_index)/2 );
        
        % ---- Check if iteration limit is reached ----
        counter = counter + 1;
        if counter >= iteration_limit
            err_msg = 'No root found after iteration limit';
            error(err_msg)
        end
    end
end