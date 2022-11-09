function [online_parameters] = gen_proposal(offline_data)
    % online_parameters_multi = { G, mu, rho, rc, [mass0, mass1, mass2]}
    % offline_data.online_parameters_property{1} - range
    % offline_data.online_parameters_property{2} - mean
%     G       = offline_data.online_parameters_property{2}{1} + (rand( size(offline_data.online_parameters_property{2}{1}) ) - 0.5).*offline_data.online_parameters_property{1}{1};
%     mu      = offline_data.online_parameters_property{2}{2} + (rand( size(offline_data.online_parameters_property{2}{2}) ) - 0.5).*offline_data.online_parameters_property{1}{2};
%     rho     = offline_data.online_parameters_property{2}{3} + (rand( size(offline_data.online_parameters_property{2}{3}) ) - 0.5).*offline_data.online_parameters_property{1}{3};
%     rc      = offline_data.online_parameters_property{2}{4} + (rand( size(offline_data.online_parameters_property{2}{4}) ) - 0.5).*offline_data.online_parameters_property{1}{4};
%     mass    = offline_data.online_parameters_property{2}{5} + (rand( size(offline_data.online_parameters_property{2}{5}) ) - 0.5).*offline_data.online_parameters_property{1}{5};
%     online_parameters_multi = {G, mu, rho, rc, mass};


% ------ online_parameters = [ G, mu, rho, rc, mass0 ] ------
    % offline_data.online_parameters_property{1} - range
    % offline_data.online_parameters_property{2} - mean
    % online_parameters = offline_data.online_parameters_property{2} + (rand(1,5) - 0.5).*offline_data.online_parameters_property{1};

%------- online_parameters = [ G, mu, rho, rc, mass0 ] ---------
    % G, mu, rho, rc - uniformly distributed 
    % mass0 - log-uniformly distribution in 10^(4.2) to 10^(7)
    % offline_data.online_parameters_property{1} - range
    % offline_data.online_parameters_property{2} - mean
    %     log_mass0_mean = (4.2+7)/2;
    %     log_mass0_range = 7-4.2;
    %     prior_mean = [G, mu, rho, rc, log_mass0_mean];
    %     op_range = [ 0.2*abs([G, mu, rho, rc]) , log_mass0_range ];
    %     op_mean = prior_mean;
    %     offline_data.online_parameters_property = { op_range, op_mean };

    % Generate for G, mu, rho and rc first
    online_parameters = offline_data.online_parameters_property{2}(1:4) + (rand(1,4) - 0.5).*offline_data.online_parameters_property{1}(1:4);
    
    % Generate the log-uniform distribution of mass0
    online_parameters(5) = 10^( offline_data.online_parameters_property{2}(5) + (rand-0.5)*offline_data.online_parameters_property{1}(5) );

end