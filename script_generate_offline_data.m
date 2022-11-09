load('Online_variables')
G = online_data(1);
mu = online_data(2);
rho = online_data(3);
rc = online_data(4);

% ---- Initializing offline data.constants ------
    P_atm= 1e+5; %atmospheric pressure
    g= 9.81; %gravitational acceleration
    nu = 0.25; %Poisson's ratio
    dist = 500; %radial distance from the top of the source
    Rgas= 462; % J*kg-1*K-1
    Temperature = 1100; % in Kelvin K = C+273.15  for magma = 1100 for Lab analogue = 289.1500
    InitialMagmaLevel =350; %depth magma free surface relative to groud surface Stromboli = 50m
    Ls0 = nan; % Bubble initial length Stromboli = 5m


% ----- Initialize variables for quadrature -------------------
    order_of_ODE=2;
    % no_of_RKGL_nodes=1;

    [ fun_list ] = fun_RHS_f_explicitRK( order_of_ODE );
    % [ dfun_list ] = dfun_RHS_f_dLj( order_of_ODE );
    % [ RK_A, RK_b, RK_c ] = Runge_Kutta_butchertableau( no_of_RKGL_nodes );

        offline_data.quadrature = {fun_list};
        offline_data.constants = {P_atm, g, nu, dist, Rgas, Temperature, InitialMagmaLevel, Ls0, order_of_ODE};


% ------- Defining the mean and range of prior ---------
    % G, mu, rho, rc - uniformly distributed 
    % mass0 - log-uniformly distribution in 10^(4.2) to 10^(7)
    log_mass0_lower = 5.2;
    log_mass0_upper = 7;
    log_mass0_mean = (log_mass0_lower+log_mass0_upper)/2;
    log_mass0_range = log_mass0_upper-log_mass0_lower;
    prior_mean = [G, mu, rho, rc, log_mass0_mean];
    op_range = [ 0.2*abs([G, mu, rho, rc]) , log_mass0_range ];
    op_mean = prior_mean;

        offline_data.online_parameters_property = { op_range, op_mean };
        
save('Offline_variables','offline_data')