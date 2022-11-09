function [ E_MLMCMC, ube_pop_var_time_to_erupt, ube_pop_var_final_tilt ] = MLMCMC_Semeru_uniform_prior_explicitRK_1burn_col_timetilt( offline_data, online_parameters, delta, L, rate, starting_level, real_time_step, alpha, burn )
%   delta           - Noisy observation
%   L               - Highest level of MLMCMC
%   rate            - Convergence rate of MLMCMC O(2^{-rate*L}). Max rate is 2*no_of_RKGL_nodes
%   starting_level  - 
%   alpha           - (l+l')^{\alpha}, to reduce the L^2 factor of MLMCMC
%   beta            - For pCN
%   burn            - No. of samples to discard initially 
%   offline_data.constants = {P_atm, g, nu, dist, Rgas, Temperature, InitialMagmaLevel, Ls0, order_of_ODE, no_of_RKGL_nodes, obs_cov_mat}
%   offline_data.quadrature = {fun_list, dfun_list, RK_A, RK_b, RK_c}
%   offline_data.online_parameters_property = { op_range, op_mean, op_percent_range_to_move }
%   online_parameters = [G, mu, rho, mass0, rc]


% ----- Initialize variables for quadrature -------------------
% order_of_ODE = offline_data.constants{9};
% 
% [ fun_list ] = fun_RHS_f_explicitRK( order_of_ODE );
% % [ dfun_list ] = dfun_RHS_f_dLj( order_of_ODE );
% % [ RK_A, RK_b, RK_c ] = Runge_Kutta_butchertableau( no_of_RKGL_nodes );
% % offline_data.quadrature = {fun_list, dfun_list, RK_A, RK_b, RK_c};
% offline_data.quadrature = {fun_list};
% ------------------------------------------------------

% ----- Initial Guess of Parameters to be found using MCMC -----
% G= 10^11;%shear modulus Stromboli = 10^9.4; Gelatin = 5000
% mu= 10^(3.4); %viscosity Stromboli 10^3
% rho= 2400; %magma density Stromboli = 2600
% rc =25; %conduit radius Stromboli = 3m
% mass0 = 10^6;%5,5.5,5.7,6; %kg Stromboli = 780 for analogue experiment %symbol "n" in the poster

% online_parameters = [G, mu, rho, rc, mass0];
% ----------------------------------------------------

% ----- Compute the eigenvectors and eigenvalues for covariance matrix ----
% [obs_cov_evectors, obs_cov_evalues] = eig(obs_cov_mat);
% prior_cov_mat = offline_data.constants{12};
% [prior_cov_evectors, prior_cov_evalues] = eig(prior_cov_mat);
% offline_data.constants{end+1} = prior_cov_evectors;
% offline_data.constants{end+1} = prior_cov_evalues;


% ---- 12/02/2020 --------
col_time_to_erupt = [];
col_final_tilt = [];
% -------------------

if L < starting_level
    disp('We are starting at L=',num2str(starting_level))
    return;
else
    E_MLMCMC = 0;
    
    MCMC_time_level = starting_level;

    [ initial_old_online_parameters_wrt_l ] = Burn_in_lognormal( offline_data, online_parameters, real_time_step, MCMC_time_level, delta, burn );
    
    % ***** When l' = starting_level ******
%     l_prime = starting_level;
    
    %[l, l_prime] %#######################
    
    % Expectation of B^{1}(u) wrt \rho^{J_{1}, 1, \delta}
    [ E_B, col_time_to_erupt_part, col_final_tilt_part ] = Plain_MCMC_B_wrt_rho_l_when_l_and_lprime_are_1( offline_data, real_time_step, MCMC_time_level, delta, L, rate, alpha, initial_old_online_parameters_wrt_l );
    E_MLMCMC = E_MLMCMC + E_B;
    col_time_to_erupt = [col_time_to_erupt; col_time_to_erupt_part];
    col_final_tilt = [col_final_tilt; col_final_tilt_part];
    % ************************
    

%     % ------***** When l'>=starting_level+1 *******------
%     for l_prime = starting_level+1:L
%         %[l, l_prime] %#######################
%         
%         
%         % Expectation of B^{l'}(u) wrt \rho^{J_{1}, 1, \delta}, for l' = 2:L-l
%         E_B = Plain_MCMC_B_wrt_rho_l_when_l_is_1( matrix_F, vector_Obs, vector_QoI, osci_c, L_nodes, L_weights, left, right, delta, L, MCMC_time_level, l_prime, rate, alpha, beta, initial_old_u_wrt_l );
%         E_MLMCMC = E_MLMCMC + E_B;
%     end
%     % ------***********************------
    
    initial_old_online_parameters_wrt_ln1 = initial_old_online_parameters_wrt_l;
    
    % --------------- When l >= starting_level+1 ----------------------------------------------------
    for MCMC_time_level = starting_level+1:L       % the upper limit should be able to leave it at L
        
%         [ initial_old_online_parameters_wrt_l ] = Burn_in_lognormal( offline_data, online_parameters, real_time_step, MCMC_time_level, delta, burn );
        
%         l_prime = starting_level;
        
        %Plain_MCMC_rho_l_when_lprime_is_1 Computes the estimated expectation of
        %A_{1}^{l}(u)*B^{1}(u)*I^{l}(u), A_{1}^{l}(u)*I^{l}(u), 
        %(1 - A_{1}^{l}(u))*B^{1}(u)*I^{l}(u), B^{1}(u)*( 1-I^{l}(u) )
        %with respect to the posterior probability measure \rho^{J_l, l, \delta}
        [ E_l2_A1_B_I, E_l2_A1_I, E_l2_1nA1_B_I, E_l2_B_1nI, col_time_to_erupt_part, col_final_tilt_part ] = Plain_MCMC_rho_l_when_lprime_is_1( offline_data, real_time_step, MCMC_time_level, delta, L, rate, alpha, initial_old_online_parameters_wrt_l );
        col_time_to_erupt = [col_time_to_erupt; col_time_to_erupt_part];
        col_final_tilt = [col_final_tilt; col_final_tilt_part];
        [ E_l1_A2_B_1nI, E_l1_A2_1nI, E_l1_1nA2_B_1nI, E_l1_B_I, col_time_to_erupt_part, col_final_tilt_part ] = Plain_MCMC_rho_lminus1_when_lprime_is_1( offline_data, real_time_step, MCMC_time_level, delta, L, rate, alpha, initial_old_online_parameters_wrt_ln1 );
        col_time_to_erupt = [col_time_to_erupt; col_time_to_erupt_part];
        col_final_tilt = [col_final_tilt; col_final_tilt_part];
        E_MLMCMC = E_MLMCMC + ( E_l2_A1_B_I - E_l2_A1_I*E_l1_B_I + E_l1_A2_1nI*E_l2_1nA1_B_I ) - ( E_l1_A2_B_1nI - E_l1_A2_1nI*E_l2_B_1nI + E_l2_A1_I*E_l1_1nA2_B_1nI);
        
%         % Expectation of A^{l}(u)*B^{1}(u) <E_AB> and -A^{l}(u) <E_nA> wrt \rho^{J_{l}, l, \delta}
%         [ E_AB, E_nA ] = Plain_MCMC_AB_nA_when_lprime_is_1( matrix_F, vector_Obs, vector_QoI, osci_c, L_nodes,L_weights,left, right, delta, L, l, l_prime, rate, beta, initial_old_u_wrt_l );
%         
%         % Expectation of B^{1}(u) wrt \rho^{J_{l-1}, l-1, \delta}
%         E_B = Plain_MCMC_B_wrt_rho_ln1_when_lprime_is_1( matrix_F, vector_Obs, vector_QoI, osci_c, L_nodes,L_weights,left, right, delta, L, l, l_prime, rate, beta, initial_old_u_wrt_ln1 );
%         E_MLMCMC = E_MLMCMC + E_AB + E_nA*E_B;
        
%         for l_prime = starting_level+1:L-MCMC_time_level
%             
%             [ E_l2_A1_B_I, E_l2_A1_I, E_l2_1nA1_B_I, E_l2_B_1nI ] = Plain_MCMC_rho_l( matrix_F, vector_Obs, vector_QoI, osci_c, L_nodes, L_weights, left, right, delta, L, MCMC_time_level, l_prime, rate, alpha, beta, initial_old_u_wrt_l );
%             [ E_l1_A2_B_1nI, E_l1_A2_1nI, E_l1_1nA2_B_1nI, E_l1_B_I ] = Plain_MCMC_rho_lminus1( matrix_F, vector_Obs, vector_QoI, osci_c, L_nodes, L_weights, left, right, delta, L, MCMC_time_level, l_prime, rate, alpha, beta, initial_old_online_parameters_wrt_ln1 );
%             E_MLMCMC = E_MLMCMC + ( E_l2_A1_B_I - E_l2_A1_I*E_l1_B_I + E_l1_A2_1nI*E_l2_1nA1_B_I ) - ( E_l1_A2_B_1nI - E_l1_A2_1nI*E_l2_B_1nI + E_l2_A1_I*E_l1_1nA2_B_1nI);
%             
% %             % Expectation of A^{l}(u)*B^{l'}(u) <E_AB> and -A^{l}(u) <E_nA> wrt \rho^{J_{l}, l, \delta}
% %             [ E_AB, E_nA ] = Plain_MCMC_AB_nA( matrix_F, vector_Obs, vector_QoI, osci_c, L_nodes,L_weights,left, right, delta, L, l, l_prime, rate, beta, initial_old_u_wrt_l );
% %             
% %             % Expectation of B^{l'}(u) wrt \rho^{J_{l-1}, l-1, \delta}
% %             E_B = Plain_MCMC_B_wrt_rho_ln1( matrix_F, vector_Obs, vector_QoI, osci_c, L_nodes,L_weights,left, right, delta, L, l, l_prime, rate, beta, initial_old_u_wrt_ln1 );
% %             E_MLMCMC = E_MLMCMC + E_AB + E_nA*E_B;
%         end
        initial_old_online_parameters_wrt_ln1 = initial_old_online_parameters_wrt_l;
    end
end

mean_time_to_erupt = mean(col_time_to_erupt);
ube_pop_var_time_to_erupt =  sum((col_time_to_erupt-mean_time_to_erupt).^2)/(length(col_time_to_erupt)-1);
mean_final_tilt = mean(col_final_tilt);
ube_pop_var_final_tilt = sum((col_final_tilt-mean_final_tilt).^2)/(length(col_final_tilt)-1);

E_MLMCMC = E_MLMCMC';
E_MLMCMC = E_MLMCMC(:)';
end



%% --- Edit this for different ODE (part 1) -----------------------
%
%
%
function [ fun_list ] = fun_RHS_f_explicitRK( order_of_ODE )
%fun_RHS_f Right hand side of a ODE. If order of ODE is >1, fun_RHS_f is
%the vector form, of length equal to order, after changing problem into 1st order ODE.
%   order - order of ODE
%   return a cell array of functions as fun_list

fun_list = cell(1,order_of_ODE);
var_names = ['(d1, d2, d3, d4, d_h1, u_s, A, t, L)'];
% for i = 1:order
%     new_var = [',L',num2str(i)];
%     var_names = [var_names, new_var];
% end
% var_names = [var_names, ')'];

for i = 1:order_of_ODE
    function_name = ['fun',num2str(i)];
    init_syntax = ['fun_list{',num2str(i),'}=@', var_names, function_name, var_names, ';'];
    eval(init_syntax)
end

end

function value1 = fun2(d1, d2, d3, d4, d_h1, u_s, A, t, L)
% value1 = d1*(  1./(L{1}.*(d_h1 - A*L{1} - u_s.*t)) ) + d2 + d3*(  1./(d_h1 - A*L{1} - u_s.*t)  ) + d4*L{2} + fun_R(d1, d2, d3, d4, d_h1, u_s, A, t);
value1 = d1*(  1./(L(1).*(d_h1 - A*L(1) - u_s.*t)) ) + d2 + d3*(  1./(d_h1 - A*L(1) - u_s.*t)  ) + d4*L(2);
% value1 = L{2} + 6*L{1};
end

% function value = fun_R(d1, d2, d3, d4, d_h1, u_s, A, t)
% C = pi/10;
% L = @(t) exp(-cos(C*t));
% value = (C^2)*L(t).*(cos(C*t) + (sin(C*t)).^2) - d1*(  1./(L(t).*(d_h1 - A*L(t) - u_s.*t)) ) - d2 - d3*(  1./(d_h1 - A*L(t) - u_s.*t)  ) - d4*C*sin(C*t).*L(t);
% end

function value2 = fun1(d1, d2, d3, d4, d_h1, u_s, A, t, L)
value2 = L(2);
end
%
%
%
% -----------------------------------------------------------------

%% --- Edit this for different ODE (part 2) -----------------------
%
%
%
function [ dfun_list ] = dfun_RHS_f_dLj( order_of_ODE )
%dfun_RHS_f_dLj all partial derivatives of Right hand side of a ODE with respect to different Lj. 
%   order - order of ODE
%   return a cell array of functions as fun_list

dfun_list = cell(order_of_ODE,order_of_ODE);
var_names = '(d1, d2, d3, d4, d_h1, u_s, A, t, L)';
% for i = 1:order
%     new_var = [',L',num2str(i)];
%     var_names = [var_names, new_var];
% end
% var_names = [var_names, ')'];

for i = 1:order_of_ODE
    for j = 1:order_of_ODE
        function_name = ['dfun',num2str(i),num2str(j)];
        init_syntax = ['dfun_list{',num2str(i),',',num2str(j),'}=@', var_names, function_name, var_names, ';'];
        eval(init_syntax)
    end
end

end

function value1 = dfun21(d1, d2, d3, d4, d_h1, u_s, A, t, L)
value1 = -d1*((d_h1*L{1} - A*L{1}.^2 - u_s*t.*L{1}).^(-2)).*(d_h1 - 2*A*L{1} - u_s.*t) + d3*( (d_h1 - A*L{1} - u_s.*t).^(-2) )*A;
% value1 = 6*ones( length(L{1}), 1);
end

function value1 = dfun22(d1, d2, d3, d4, d_h1, u_s, A, t, L)
value1 = d4*ones( length(L{1}) ,1);
% value1 = 1*ones( length(L{1}) ,1);
end

function value2 = dfun11(d1, d2, d3, d4, d_h1, u_s, A, t, L)
value2 = 0*ones( length(L{1}) ,1);
end

function value2 = dfun12(d1, d2, d3, d4, d_h1, u_s, A, t, L)
value2 = 1*ones( length(L{1}) ,1);
end
%
%
%
% -----------------------------------------------------------------

%% --- Function to generate RK butcher tableau -----
function [ RK_A, RK_b, RK_c ] = Runge_Kutta_butchertableau( no_of_RKGL_nodes )
%UNTITLED6 Summary of this function goes here
%   s - no. of GL points in each time interval

[vector_GL_points,w]=gauleg(no_of_RKGL_nodes);
RK_A = zeros(no_of_RKGL_nodes,no_of_RKGL_nodes);
degree_poly = no_of_RKGL_nodes-1;
a = -1;
for i = 1:length(vector_GL_points)
    RK_b = vector_GL_points(i);
    new_nodes = vector_GL_points*(RK_b-a)/2 + (a+RK_b)/2;
    for j = 0:degree_poly
        RK_A(i,j+1) = ((RK_b-a)/2)*Legendre_polynomial( vector_GL_points, j, new_nodes )*w';
    end
end
RK_A = RK_A/2;

% b is a row vector
RK_b = w/2;

RK_c = (vector_GL_points+1)/2;

end

function [ value_at_x ] = Legendre_polynomial( vector_GL_points, i, x )
%Legendre_polynomial Generate the Legendre polynomial at the point x, with
%degree = length(vector_GL_points), in the interval [-1,1]
%   vector_GL_points - vector containing the points of the Legendre polynomial in the interval [-1,1]
%   i - ith polynomial (goes from 0 to n, the degree)
%   x - point to be evaluated

degree_poly = length(vector_GL_points)-1;
value_at_x = 1;
for k = 0:(i-1)
    value_at_x = value_at_x .* (x - vector_GL_points(k+1))./(vector_GL_points(i+1)-vector_GL_points(k+1));     %we use k+1 because matlab indicing starts from 1
end
for k = (i+1):degree_poly
    value_at_x = value_at_x .* (x - vector_GL_points(k+1)) ./(vector_GL_points(i+1)-vector_GL_points(k+1));     %we use k+1 because matlab indicing starts from 1
end

end

function [x,w]=gauleg(n)

x=zeros(n,1);
w=zeros(n,1);
m=(n+1)/2;
xm=0.0;
xl=1.0;
for i=1:m
    z=cos(pi*(i-0.25)/(n+0.5));
    while 1
        p1=1.0;
        p2=0.0;
        for j=1:n
            p3=p2;
            p2=p1;
            p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
        end
        pp=n*(z*p1-p2)/(z*z-1.0);
        z1=z;
        z=z1-p1/pp;
        if (abs(z-z1)<eps), break, end
    end
    x(i)=xm-xl*z;
    x(n+1-i)=xm+xl*z;
    w(i)=2.0*xl/((1.0-z*z)*pp*pp);
    w(n+1-i)=w(i);
end
x=x';
w=w';
return;
end
% --------------------------------------------------------

%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$ --- Function to compute observation ----- $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
function [ tilt_total, time_of_erupt ] = MLMCMC_fun_observation_explicit_RK( offline_data, real_time_step, MCMC_time_level, online_parameters )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%   offline_data.constants = {P_atm, g, nu, dist, Rgas, Temperature, InitialMagmaLevel, Ls0, order_of_ODE, no_of_RKGL_nodes, obs_cov_mat, prior_cov_mat, prior_cov_evectors, prior_cov_evalues}
%   offline_data.quadrature = {fun_list, dfun_list, RK_A, RK_b, RK_c}
%   online_parameters = {G, mu, rho, mass0, rc}

MCMC_time_step = 4*real_time_step*(2^(-MCMC_time_level));

% % ----- Parameters to be found using MCMC -----
    % G= 10^11;%shear modulus Stromboli = 10^9.4; Gelatin = 5000
    % mu= 10^3.4; %viscosity Stromboli 10^3
    % rho= 2400; %magma density Stromboli = 2600
    % mass0 = 10^6;%5,5.5,5.7,6; %kg Stromboli = 780 for analogue experiment %symbol "n" in the poster
    % rc =25; %conduit radius Stromboli = 3m
    G = online_parameters(1);
    mu = online_parameters(2);
    rho = online_parameters(3);
    rc = online_parameters(4);
    mass0 = online_parameters(5);
% % ----------------------------------------------------

% % ----- Constants -----
    % P_atm= 1e+5; %atmospheric pressure
    % g= 9.81; %gravitational acceleration
    % nu = 0.25; %Poisson's ratio
    % dist = 500; %radial distance from the top of the source
    % 
    % % PV = nRT
    % Rgas= 462; % J*kg-1*K-1
    % Temperature = 1100; % in Kelvin K = C+273.15  for magma = 1100 for Lab analogue = 289.1500
    % 
    % % time_frame =0.5; %secs
    % InitialMagmaLevel =350; %depth magma free surface relative to groud surface Stromboli = 50m 
    % Ls0 = nan; % Bubble initial length Stromboli = 5m
    % sigma = 0.07; %N/m surface tension ;Silicon Oil = 0.025; water = 0.07 magma = 0.07
    P_atm = offline_data.constants{1};
    g = offline_data.constants{2};
    nu = offline_data.constants{3};
    dist = offline_data.constants{4};
    Rgas = offline_data.constants{5};
    Temperature = offline_data.constants{6};
    InitialMagmaLevel = offline_data.constants{7};
    Ls0 = offline_data.constants{8};
    order_of_ODE = offline_data.constants{9};
% % -----------------------

% % -- Quadrature constants -----------
    % offline_data.quadrature = {fun_list, dfun_list, RK_A, RK_b, RK_c}
    fun_list = offline_data.quadrature{1};
%     dfun_list = offline_data.quadrature{2};
%     RK_A = offline_data.quadrature{3};
%     RK_b = offline_data.quadrature{4};
%     RK_c = offline_data.quadrature{5};
% % ------------------------------------

% ----- (Dependent Constants) Initial condition for non-time-independent variables ----
    Nf= (rho/mu)*sqrt(g*8*(rc^3)); %buoyancy Reynolds number (dimensionless inverse viscosity)
    Fr= 0.34*((1+((31.08/Nf)^1.45))^(-0.71)); %Froude number (dimensionless velocity) that rapresent the ratio of inertial and gravitational forces
    lbd= (6*(Fr/Nf))^(1/3); %dimensionless falling film thickness  (lambda' in the poster)
    lbd1= lbd*rc; %falling film thickness  (lambda in the poster)
    rs=rc-lbd1; %gas slug radius
    Vb = (2*(rho)*g*(lbd1^3))/(3*mu*rs); % bubble raising velocity m/s  (u_b in poster)
% ----------------------------------------------------------------

% ----- (Dependent Constants) Initial condition for time-dependent variables -----
    % Computed from the Constants above
    if isnan(Ls0)
        Ls0 =rs*2;%gas slug initial length
    end
    VolInitial= (pi*(rs^2)*Ls0); %Initial bubble volume in m^3
    PressureInitial = (mass0*Rgas*Temperature)/VolInitial;  %mass0 is "n" in the poster, PressureInitial is NOT P_0 in poster
    DepthInitial = (PressureInitial-P_atm)/(rho*g); %bubble head depth. DepthInitial is h_0 in poster
    % I think DepthInitial should be (PressureInitial-P_atm)/(rho*g)
    % disp('DepthInitial might not be correct')
% ---------------------------------------------------------------


% ----- (Dependent Constants) Initializing constants used in the ODE -----
    A1 = 1-(rs^2/rc^2); %dimensionless parameter related with the conduid cross-sectional area occupied by the magma     "A' in the poster"

% mass_limit = (4*A1*pi*(rs^4)*rho*g)/(Rgas*Temperature)



% ----- (Compute height limit of magma) -----
    La= (PressureInitial*Ls0)/P_atm;  %length that the slug would have at atmospheric pressure
    Ps_lim= sqrt(rho*g*A1*P_atm*La); %pressure after that the slug became instable from DelBello (2011), HAVE THE SAME VALUE OF JAMES (2009)
    h_min_height_of_magma_before_erupt = ((Ps_lim -P_atm)/(rho*g)); 

% ---- Compute Altitude of station -----
    StationAltitude = DepthInitial+Ls0+InitialMagmaLevel; %Bubble Base deph with reference to the magma top level  %%%%%%%%%%%%%%%%%%%%StationAltitude
    % RimAltitude = DepthInitial+Ls0+abs(InitialMagmaLevel); %Bubble Base deph with reference to the magma top level  %%%%%%%%%%%%%%%%%%%%StationAltitude
    
    
% ---- Computing constants used in ODE ----
    P0 = PressureInitial; %DepthInitial*rho*g + P_atm;
    d1 = 2/(rho*(2-A1))*P0*Ls0;
    d2 = -2*g/(2-A1);
    d3 = -2*P_atm/(rho*(2-A1));
    d4 = -16*mu/(rho*(2-A1)*(rc^2));

    d_h1 = DepthInitial + Ls0*A1;
% ----------------------------------------------------------


% ------ Initializing initial condition for ODE -----------------
    t_0 = 0;
    L_0 = [Ls0; 0];
% --------------------------------------------------------------


% ---- Max iteration needed base on the time the bottom of the bubble reach
    t_lim =(StationAltitude)/Vb; %time the slug reach the surface                   %%%%%%%%%%%%%%%%%%%%StationAltitude
    iteration_limit = ceil(t_lim/MCMC_time_step);
% --------------------------------------------------



% ----- Solving for the bubble length at custom time nodes ---------------
% tic
[ L_collect, h_collect, s_collect, t_collect ] = MLMCMC_RK_explicit_bubble_length_fix_step( order_of_ODE, fun_list, t_0, MCMC_time_step, iteration_limit, L_0, d1, d2, d3, d4, d_h1, Vb, A1, h_min_height_of_magma_before_erupt );
% 'bubble'
% toc
% -----------------------------------------------------------------
% tic
[interpolated_L_collect, interpolated_dL_dt_collect, interpolated_h_collect, interpolated_s_collect, interpolated_t_collect] = interpolate_time_from_h( L_collect, h_collect, t_collect, Vb, A1, d1, d2, d3, d4, d_h1, h_min_height_of_magma_before_erupt );
time_of_erupt = interpolated_t_collect(end);
% 'interpolate'
% toc
% ----- Position of (z1)Top of magma, (z2)Top of bubble, (z3)Bottom of bubble, (z4)Conduit bed -------------------------
% z4 = StationAltitude;           % Depth of chamber wrt the top of the conduit
% z3 = z4 - s_collect;            % Position of the bubble base 
% z2 = z3 - L_collect(:,1);       % Position of the bubble top
% z1 = z2 - h_collect;            % Position of the magma top
% if abs(InitialMagmaLevel - z1(1))<1e-12
%     z0 = InitialMagmaLevel;     % initial magma position
% else
%     InitialMagmaLevel
%     z1(1)
%     InitialMagmaLevel - z1(1)
%     disp('Something is wrong, Initial magma level should be equal to z1^{(0)}');
%     return
% end

z4 = StationAltitude;           % Depth of chamber wrt the top of the conduit
z3 = z4 - interpolated_s_collect;            % Position of the bubble base 
z2 = z3 - interpolated_L_collect;       % Position of the bubble top
z1 = z2 - interpolated_h_collect;            % Position of the magma top
% z4 = StationAltitude;           % Depth of chamber wrt the top of the conduit
% z3 = z4 - s_collect;            % Position of the bubble base 
% z2 = z3 - L_collect(:,1);       % Position of the bubble top
% z1 = z2 - h_collect;            % Position of the magma top
z0 = InitialMagmaLevel;     % initial magma position

% ------------------------------------------------------------------


% ---- Computing Tilt sigma, caused by normal stress ---------------
% tic
[ tilt_sigma ] = Tilt_sigma( dist, rho, g, G, nu, rc, A1, z0, z1, z2, z3, z4);
% 'tile_sigma'
% toc
% ------------------------------------------------------------------

% ----- Computing Tilt tau, caused by shear stress ---------------
% tic
[ tilt_tau ] = Tilt_tau( dist, rho, g, G, nu, mu, rc, rs, z1, z2, z3, interpolated_dL_dt_collect );
% [ tilt_tau ] = Tilt_tau( dist, rho, g, G, nu, mu, rc, rs, z1, z2, z3, L_collect(:,2) );
% 'tilt_tau'
% toc
% ----------------------------------------------------------------

tilt_total = tilt_sigma + tilt_tau;

end

    % ##################################################################################################################################################
    function [ tilt_sigma ] = Tilt_sigma( dist, rho, g, G, nu, rc, A1, z0, z1, z2, z3, z4)

        tilt_sigma = zeros(length(z1),1);
        for i = 1:length(z1)
            tilt_sigma_at_t = 0;

            % We let \Delta P = mz + c
            if z0 < z2(i)
                % In [z1(i), z0]
                delta_P_m = rho*g;
                delta_P_c = -rho*g*z1(i);
                tilt_sigma_at_t = tilt_sigma_at_t + fun_tilt_sigma( rc, G, nu, z1(i), z0, dist, delta_P_m, delta_P_c );


                % In [z0, z2(i)]
                delta_P_m = 0;
                delta_P_c = rho*g*(z0 - z1(i));
                tilt_sigma_at_t = tilt_sigma_at_t + fun_tilt_sigma( rc, G, nu, z0, z2(i), dist, delta_P_m, delta_P_c );

                % In [z2(i), z3(i)]
                delta_P_m = rho*g*(A1-1);
                delta_P_c = rho*g*( z0 - z1(i) - (A1-1)*z2(i) );
                tilt_sigma_at_t = tilt_sigma_at_t + fun_tilt_sigma( rc, G, nu, z2(i), z3(i), dist, delta_P_m, delta_P_c );

                % In [z3(i), z4]
                delta_P_m = 0;
                delta_P_c = rho*g*( (A1-1)*z3(i) + z0 - z1(i) - (A1-1)*z2(i) );
                tilt_sigma_at_t = tilt_sigma_at_t + fun_tilt_sigma( rc, G, nu, z3(i), z4, dist, delta_P_m, delta_P_c );


            elseif (z2(i)<z0) && (z0<z3(i))
    %             disp(2);
                % In [z1(i), z2(i)]
                delta_P_m = rho*g;
                delta_P_c = -rho*g*z1(i);
                tilt_sigma_at_t = tilt_sigma_at_t + fun_tilt_sigma( rc, G, nu, z1(i), z2(i), dist, delta_P_m, delta_P_c );

                % In [z2(i), z0]
                delta_P_m = rho*g*A1;
                delta_P_c = rho*g*( (1-A1)*z2(i) - z1(i));
                tilt_sigma_at_t = tilt_sigma_at_t + fun_tilt_sigma( rc, G, nu, z2(i), z0, dist, delta_P_m, delta_P_c );


                % In [z0, z3(i)]
                delta_P_m = rho*g*(A1-1);
                delta_P_c = rho*g*( z0 - z1(i) + (1-A1)*z2(i) );
                tilt_sigma_at_t = tilt_sigma_at_t + fun_tilt_sigma( rc, G, nu, z0, z3(i), dist, delta_P_m, delta_P_c );

                % In [z3(i), z4]
                delta_P_m = 0;
                delta_P_c = rho*g*( (A1-1)*z3(i) + z0 - z1(i) + (1-A1)*z2(i) );
                tilt_sigma_at_t = tilt_sigma_at_t + fun_tilt_sigma( rc, G, nu, z3(i), z4, dist, delta_P_m, delta_P_c );

            elseif (z3(i)<z0) && (z0<z4)
    %             disp(3);
                % In [z1(i), z2(i)]
                delta_P_m = rho*g;
                delta_P_c = -rho*g*z1(i);
                tilt_sigma_at_t = tilt_sigma_at_t + fun_tilt_sigma( rc, G, nu, z1(i), z2(i), dist, delta_P_m, delta_P_c );

                % In [z2(i), z3(i)]
                delta_P_m = rho*g*A1;
                delta_P_c = rho*g*( (1-A1)*z2(i) - z1(i));
                tilt_sigma_at_t = tilt_sigma_at_t + fun_tilt_sigma( rc, G, nu, z2(i), z3(i), dist, delta_P_m, delta_P_c );


                % In [z3(i), z0]
                delta_P_m = rho*g;
                delta_P_c = rho*g*( (1-A1)*z2(i) + (A1-1)*z3(i) - z1(i) );
                tilt_sigma_at_t = tilt_sigma_at_t + fun_tilt_sigma( rc, G, nu, z3(i), z0, dist, delta_P_m, delta_P_c );

                % In [z0, z4]
                delta_P_m = 0;
                delta_P_c = rho*g*( (A1-1)*z3(i) + z0 - z1(i) + (1-A1)*z2(i) );
                tilt_sigma_at_t = tilt_sigma_at_t + fun_tilt_sigma( rc, G, nu, z0, z4, dist, delta_P_m, delta_P_c );

            else
                disp('Position of initial magma level cannot be greater than position of conduit bed')
                return
            end
            tilt_sigma(i) = tilt_sigma_at_t;
        end

    end
    % -------------------------------------------------------
    function [ tilt ] = fun_tilt_sigma( rc, G, mu, lower_limit, upper_limit, station_radius_dist, delta_P_m, delta_P_c )
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here

    % tilt_sigma = -dU_{\sigma}/dr
    tilt = -((rc^2)/G)*( delta_P_m*( df_1_mz( mu, lower_limit, upper_limit, station_radius_dist ) - 1.5*df_2_mz( lower_limit, upper_limit, station_radius_dist ) ) + delta_P_c*( df_1_c( lower_limit, upper_limit, station_radius_dist ) + df_2_c( mu, lower_limit, upper_limit, station_radius_dist ) ) );

    end
    % -------------------------------------------------------
    function [ value ] = df_1_mz( nu, lower_limit, upper_limit, station_radius_dist )
    %df_1_mz one part of the derivative of the normal stress 
    %   nu - Poisson's ratio
    %   lower_limit - lower limit of the integral, which depends on time
    %   upper_limit - upper limit of the integral, which depends on time
    %   station_radius_dist - the distance of the station from the center of the conduit at the
    %   surface

    %     fun_d_f_1_mz = @(r,z) (nu + 1)*( (1./(z+sqrt(r.^2 + z.^2))).*(r./(sqrt(r.^2 + z.^2))) - 1./r + (r.*z)./((r.^2 + z.^2).^(1.5)) );
    %     value = fun_d_f_1_mz(station_radius_dist, upper_limit) - fun_d_f_1_mz(station_radius_dist, lower_limit);
    r = station_radius_dist;
    z = upper_limit;
    value_upper = (nu + 1)*( (1./(z+sqrt(r.^2 + z.^2))).*(r./(sqrt(r.^2 + z.^2))) - 1./r + (r.*z)./((r.^2 + z.^2).^(1.5)) );
    z = lower_limit;
    value_lower = (nu + 1)*( (1./(z+sqrt(r.^2 + z.^2))).*(r./(sqrt(r.^2 + z.^2))) - 1./r + (r.*z)./((r.^2 + z.^2).^(1.5)) );
    value = value_upper - value_lower;


    end
    % -------------------------------------------------------
    function [ value ] = df_2_mz( lower_limit, upper_limit, station_radius_dist )
    %df_2_mz one part of the derivative of the normal stress 
    %   nu - Poisson's ratio
    %   lower_limit - lower limit of the integral, which depends on time
    %   upper_limit - upper limit of the integral, which depends on time
    %   station_radius_dist - the distance of the station from the center of the conduit at the
    %   surface

    %     fun_d_f_2_mz = @(r,z) ( (1./(z+sqrt(r.^2 + z.^2))).*(r./(sqrt(r.^2 + z.^2))) - 1./r + (r.*z.*(r.^2 + 2*z.^2))./((r.^2 + z.^2).^(2.5)) );
    %     value = fun_d_f_2_mz(station_radius_dist, upper_limit) - fun_d_f_2_mz(station_radius_dist, lower_limit);

    r = station_radius_dist;
    z = upper_limit;
    value_upper = ( (1./(z+sqrt(r.^2 + z.^2))).*(r./(sqrt(r.^2 + z.^2))) - 1./r + (r.*z.*(r.^2 + 2*z.^2))./((r.^2 + z.^2).^(2.5)) );
    z = lower_limit;
    value_lower = ( (1./(z+sqrt(r.^2 + z.^2))).*(r./(sqrt(r.^2 + z.^2))) - 1./r + (r.*z.*(r.^2 + 2*z.^2))./((r.^2 + z.^2).^(2.5)) );
    value = value_upper - value_lower;

    end
    % -------------------------------------------------------
    function [ value ] = df_1_c( lower_limit, upper_limit, station_radius_dist )
    %df_1_c one part of the derivative of the normal stress 
    %   lower_limit - lower limit of the integral, which depends on time
    %   upper_limit - upper limit of the integral, which depends on time
    %   station_radius_dist - the distance of the station from the center of the conduit at the
    %   surface

    %     fun_d_f_1_c = @(r,z) -1.5*(r.*(z.^2))./((r.^2 + z.^2).^(2.5));
    %     value = fun_d_f_1_c(station_radius_dist, upper_limit) - fun_d_f_1_c(station_radius_dist, lower_limit);

    r = station_radius_dist;
    z = upper_limit;
    value_upper = -1.5*(r.*(z.^2))./((r.^2 + z.^2).^(2.5));
    z = lower_limit;
    value_lower = -1.5*(r.*(z.^2))./((r.^2 + z.^2).^(2.5));
    value = value_upper - value_lower;

    end
    % -------------------------------------------------------
    function [ value ] = df_2_c( nu, lower_limit, upper_limit, station_radius_dist )
    %df_2_c one part of the derivative of the normal stress 
    %   nu - Poisson's ratio
    %   lower_limit - lower limit of the integral, which depends on time
    %   upper_limit - upper limit of the integral, which depends on time
    %   station_radius_dist - the distance of the station from the center of the conduit at the
    %   surface

    %     fun_d_f_2_c = @(r,z) nu*(r)./((r.^2 + z.^2).^(1.5));
    %     value = fun_d_f_2_c(station_radius_dist, upper_limit) - fun_d_f_2_c(station_radius_dist, lower_limit);

    r = station_radius_dist;
    z = upper_limit;
    value_upper = nu*(r)./((r.^2 + z.^2).^(1.5));
    z = lower_limit;
    value_lower = nu*(r)./((r.^2 + z.^2).^(1.5));
    value = value_upper - value_lower;

    end
    % ##################################################################################################################################################
%%
function [ tilt_tau ] = Tilt_tau( dist, rho, g, G, nu, mu, rc, rs, z1, z2, z3, dL_dt )

Tbubble = -rho*g*( (rs^2)/(2*rc) - rc/2 );
Tmagma = 4*mu*dL_dt*( (rs^2)/(rc^3) );

tilt_tau = (-3/(4*pi*G))*((1-2*nu)/(1-nu))*( Tbubble.*df_tau( z2, z3, dist ) + Tmagma.*df_tau( z1, z2, dist ) );

end

% -------------------------------------------------------
function value = df_tau( lower_limit, upper_limit, station_radius_dist )
%     fun_df_tau = @(r,z) (r/3).*( 5*z.^2 + 2*r.^2 )./( (z.^2 + r.^2).^(2.5) );
%     value = fun_df_tau(station_radius_dist, upper_limit) - fun_df_tau(station_radius_dist, lower_limit);

r = station_radius_dist;
z = upper_limit;
value_upper = (r/3).*( 5*z.^2 + 2*r.^2 )./( (z.^2 + r.^2).^(2.5) );
z = lower_limit;
value_lower = (r/3).*( 5*z.^2 + 2*r.^2 )./( (z.^2 + r.^2).^(2.5) );
value = value_upper - value_lower;
end
% -----------------------------------------------------------

% ----- Hermite Interpolation ---------------------------
function [interpolated_L_collect, interpolated_dL_dt_collect, interpolated_h_collect, interpolated_s_collect, interpolated_t_collect] = interpolate_time_from_h( L_collect, h_collect, t_collect, Vb, A1, d1, d2, d3, d4, d_h1, h_min_height_of_magma_before_erupt )
% Use Hermite polynomial to h function using last 2 data, and interpolate
% the time of eruption using h_min_height_of_magma_before_erupt

tolerance = 1e-5;
iteration_limit = 1000;

% -----Use Newton's method to solve for the time of eruption given h = h_min_height_of_magma_before_erupt-----
% Rate of change of height of magma above slug at the last 2 time nodes
% dh_dt = [ -Vb - L_collect(end-1,2)*A1, -Vb - L_collect(end,2)*A1 ];
dh_dt = -Vb - L_collect(:,2)*A1;

% fun_g = @(x) x - (Hermite_poly( t_collect(end-1), t_collect(end), h_collect(end-1), h_collect(end), dh_dt(end-1), dh_dt(end), x) - h_min_height_of_magma_before_erupt)./Hermite_poly_deriv( t_collect(end-1), t_collect(end), h_collect(end-1), h_collect(end), dh_dt(end-1), dh_dt(end), x);

time_erupt = Newton_method_1D( t_collect, h_collect, dh_dt, h_min_height_of_magma_before_erupt, tolerance, iteration_limit, t_collect(end-1) );
time_diff = t_collect(end) - time_erupt;
interpolated_t_collect = t_collect(2:end)-time_diff;

% L_time_erupt = Hermite_poly(t_collect(end-1), t_collect(end), L_collect(end-1,1), L_collect(end,1), L_collect(end-1,2), L_collect(end,2), time_erupt);
interpolated_L_collect = Hermite_poly(t_collect(1:(end-1)), t_collect(2:end), L_collect(1:(end-1),1), L_collect(2:end,1), L_collect(1:(end-1),2), L_collect(2:end,2), interpolated_t_collect);

% d2L_dt2_1 = d1*(  1./(L_collect(end-1,1).*(d_h1 - A1*L_collect(end-1,1) - Vb.*t_collect(end-1))) ) + d2 + d3*(  1./(d_h1 - A1*L_collect(end-1,1) - Vb.*t_collect(end-1))  ) + d4*L_collect(end-1,2);
% d2L_dt2_2 = d1*(  1./(L_collect(end,1).*(d_h1 - A1*L_collect(end,1) - Vb.*t_collect(end))) ) + d2 + d3*(  1./(d_h1 - A1*L_collect(end,1) - Vb.*t_collect(end))  ) + d4*L_collect(end,2);
% dL_dt_time_erupt = Hermite_poly(t_collect(end-1), t_collect(end), L_collect(end-1,2), L_collect(end,2), d2L_dt2_1, d2L_dt2_2, time_erupt);
d2L_dt2 = d1*(  1./(L_collect(:,1).*(d_h1 - A1*L_collect(:,1) - Vb.*t_collect(:))) ) + d2 + d3*(  1./(d_h1 - A1*L_collect(:,1) - Vb.*t_collect(:))  ) + d4*L_collect(:,2);
interpolated_dL_dt_collect = Hermite_poly(t_collect(1:(end-1)), t_collect(2:end), L_collect(1:(end-1),2), L_collect(2:end,2), d2L_dt2(1:(end-1)), d2L_dt2(2:end), interpolated_t_collect);

interpolated_h_collect = Hermite_poly(t_collect(1:(end-1)), t_collect(2:end), h_collect(1:(end-1),1), h_collect(2:end,1), dh_dt(1:(end-1),1), dh_dt(2:end,1), interpolated_t_collect);

% s_time_erupt = Vb*time_erupt;
interpolated_s_collect = Vb*interpolated_t_collect;

end

function value = fun_f(t_collect, h_collect, dh_dt, h_min_height_of_magma_before_erupt, x)
    value = (Hermite_poly( t_collect(end-1), t_collect(end), h_collect(end-1), h_collect(end), dh_dt(end-1), dh_dt(end), x) - h_min_height_of_magma_before_erupt);
end

function value = fun_g(t_collect, h_collect, dh_dt, h_min_height_of_magma_before_erupt, x)
    value = x - (Hermite_poly( t_collect(end-1), t_collect(end), h_collect(end-1), h_collect(end), dh_dt(end-1), dh_dt(end), x) - h_min_height_of_magma_before_erupt)./Hermite_poly_deriv( t_collect(end-1), t_collect(end), h_collect(end-1), h_collect(end), dh_dt(end-1), dh_dt(end), x);
end

function value = Hermite_poly(x0, x1, f_xn, f_xnp1, df_xn, df_xnp1, x)

fun_H_10 = (1 - 2*(x-x0)./(x0-x1)).*( (x-x1)./(x0-x1) ).^2;
fun_H_11 = (1 - 2*(x-x1)./(x1-x0)).*( (x-x0)./(x1-x0) ).^2;
fun_Hhat_10 = (x-x0).*( (x-x1)./(x0-x1) ).^2;
fun_Hhat_11 = (x-x1).*( (x-x0)./(x1-x0) ).^2;

value = fun_H_10.*f_xn + fun_H_11.*f_xnp1 + fun_Hhat_10.*df_xn + fun_Hhat_11.*df_xnp1;
end

function value = Hermite_poly_deriv(x0, x1, f_xn, f_xnp1, df_xn, df_xnp1, x)

fun_dH_10 = ( 6/((x1-x0)^3) )*(x-x0).*(x-x1);
fun_dH_11 = ( 6/((x0-x1)^3) )*(x-x0).*(x-x1);
fun_dHhat_10 = ( 1/((x1-x0)^2) )*(x-x1).*(3*x - x1 - 2*x0);
fun_dHhat_11 = ( 1/((x1-x0)^2) )*(x-x0).*(3*x - x0 - 2*x1);

value = fun_dH_10.*f_xn + fun_dH_11.*f_xnp1 + fun_dHhat_10.*df_xn + fun_dHhat_11.*df_xnp1;
end

function next_guess = Newton_method_1D( t_collect, h_collect, dh_dt, h_min_height_of_magma_before_erupt, tolerance, iteration_limit, initial_guess )

% next_guess = fun(initial_guess);
next_guess = fun_g(t_collect, h_collect, dh_dt, h_min_height_of_magma_before_erupt, initial_guess);

for i = 1:iteration_limit
%     if abs(next_guess-initial_guess)<tolerance
    if abs(fun_f(t_collect, h_collect, dh_dt, h_min_height_of_magma_before_erupt, next_guess))<tolerance
        return
    else
        initial_guess = next_guess;
%         next_guess = fun(initial_guess);
        next_guess = fun_g(t_collect, h_collect, dh_dt, h_min_height_of_magma_before_erupt, initial_guess);
    end
end

err_msg = ['Method did not converge after ', num2str(iteration_limit), ' iterations'];
disp(err_msg)

next_guess = min( abs(h_collect(end-1)-h_min_height_of_magma_before_erupt), abs(h_collect(end)-h_min_height_of_magma_before_erupt) );
% err_msg = ['Method did not converge after ', num2str(iteration_limit), ' iterations'];
% error(err_msg)
end

% ------------------------------------------------------
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

% @@@@@@@@@@ --- Function to compute potential --- @@@@@@@@@@@@@@@@@@@@@@@@
function [ Phi, time_to_erupt ] = MLMCMC_fun_potential( constants, delta, full_observation, time_of_erupt )
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
    full_observation_pad_zeros = [ zeros( ceil( (no_of_obs-1)*real_time_step/MCMC_time_step), 1); full_observation ];


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

function [better_Phi, better_root] = neighbor_search_better_Phi_MgR(approx_index_of_roots, MCMC_time_step, full_observation_pad_zeros, real_time_step, no_of_obs, delta, obs_cov_mat)
% When MCMC_time_step is greater than real_time_step

    % Radius of time to search (in seconds, not number of time step)
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
    avail_jump_index_left = floor(avail_time_on_left/real_time_step);
    avail_jump_index_right = floor(avail_time_on_right/real_time_step);
    
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
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% ---------- Quantity of Interest ---------------------------
function [ b_prime ] = B_lprime_when_lprime_is_1( online_parameters, time_to_erupt, tilt_total )
%B_lprime Computes B^{1}(u) = \ell(P^{J_{1},1})
%
% *Note l_prime is always 1 here

% % Coefficient of P(FEM soln) with mesh-size [2^-1]
% coeff_of_FEM_P = FEM_P( matrix_F, L_nodes, L_weights, left, right, l_prime, osci_c, u );
% 
% % B^{1}(u) = \ell(P^{J_{1},1})
% b_prime = QoI_FEM( vector_QoI, l_prime, coeff_of_FEM_P );

b_prime = [online_parameters, time_to_erupt, time_to_erupt^2, tilt_total(end), tilt_total(end)^2];

end
% -----------------------------------------------------------

%% ################################################################################################################## 
% ----- Plain_MCMC_rho_lminus1_when_lprime_is_1 ---------------------
function [ E_l1_A2_B_1nI, E_l1_A2_1nI, E_l1_1nA2_B_1nI, E_l1_B_I, col_time_to_erupt_part, col_final_tilt_part ] = Plain_MCMC_rho_lminus1_when_lprime_is_1( offline_data, real_time_step, MCMC_time_level, delta, L, rate, alpha, initial_old_online_parameters_wrt_ln1 )
%Plain_MCMC_rho_lminus1_when_lprime_is_1 Computes the estimated expectation of
%A_{2}^{l}(u)*B^{1}(u)*( 1-I^{l}(u) ), A_{2}^{l}(u)*( 1-I^{l}(u) ), 
%(1 - A_{2}^{l}(u))*B^{1}(u)*( 1-I^{l}(u) ), B^{1}(u)*I^{l}(u)
%with respect to the posterior probability measure \rho^{J_{l-1}, l-1, \delta}
%   [L_nodes, L_weights] - nodes and weights of the Gauss Legendre quad
%   left, right - left and right endpoint of the domain
%   delta - data (Computed by 1 realization of u and \vartheta)
%   L - Levels of telescopic expansions (or the finest mesh-size 2^{-L})
%   l, l_prime - current level or term in the double summation of E_L^{MLMCMC}
%   rate - rate of convergence of Observation and Quantity of Interest (if
%          the convergence is O(h^r), then rate = r)
%
% MCMC_time_level represents l from before

%       *Note l_prime == 1 here

% Collect results cummulatively
cum_E_l1_A2_B_1nI = 0;
cum_E_l1_A2_1nI = 0;
cum_E_l1_1nA2_B_1nI = 0;
cum_E_l1_B_I = 0;

l_prime = 1;
% Plain MCMC sampling size
M_ll = ceil( c_ll_rho_lminus1_when_lprime_is_1(alpha, MCMC_time_level, l_prime)*2^( 2*rate*( L-MCMC_time_level ) ) );

col_time_to_erupt_part = zeros(M_ll,1);
col_final_tilt_part = zeros(M_ll,1);

% ------- Initializing constants -------
factor_of_real_time = 4;
MCMC_time_step = factor_of_real_time*real_time_step*(2^(-MCMC_time_level));

% Number of observation
no_of_obs = length(delta);

constants = {MCMC_time_step, real_time_step, no_of_obs, offline_data.constants{10} };
% --------------------------------------

% ------ Generate initial sample of the Markov Chain ------------------------------------------------------------
old_online_parameters = initial_old_online_parameters_wrt_ln1{1};
old_time_to_erupt = initial_old_online_parameters_wrt_ln1{2};
old_observation = initial_old_online_parameters_wrt_ln1{3};
old_Phi = initial_old_online_parameters_wrt_ln1{4};

    % ***** A_{2}^l(old_u) = 1 - exp[\Phi^{J_{l-1}, l-1}(old_u;\delta) - \Phi^{J_l, l}(old_u;\delta)] *****
    % To compute \Phi^{J_l,l}(old_u;\delta)
    % Observation G^{J_l,l}(old_u)
    [old_observation_p1, old_time_of_erupt_p1] = MLMCMC_fun_observation_explicit_RK( offline_data, real_time_step, MCMC_time_level, old_online_parameters );
    % Potential function \Phi^{J_l,l}(old_u;\delta)
    [old_Phi_p1, old_time_to_erupt_p1] = MLMCMC_fun_potential( constants, delta, old_observation_p1, old_time_of_erupt_p1 );

    % Recall I^{l}(u) = { 1 if \Phi^{J_l, l}(u;\delta) - \Phi^{J_{l-1}, l-1}(u;\delta) <= 0
    %                   { 0 otherwise
    % This acts like a decision whether to use a sample a not, which will
    % be used later
    I_l = old_Phi_p1 - old_Phi;
    
    % A_{2}^{l}(old_u) = 1 - exp[\Phi^{J_{l-1}, l-1}(old_u;\delta) - \Phi^{J_l, l}(old_u;\delta)]
    A2_l = 1 - exp( -1*I_l );
    % ******************************************************************************************************

    % ########------ B^{l'}(old_u) -----##########################
    B = B_lprime_when_lprime_is_1( old_online_parameters, old_time_to_erupt, old_observation );
    % ############################################################
% ---------------------------------------------------------------------------------------------------------------------

% -------- coefficient of v^(k) in pCN -----------------
% weight_old_u = sqrt(1-beta^2);
% -------------------------------------------------------
% no_online_para = size(initial_old_online_parameters_wrt_ln1);
i = 1;
while i <= M_ll
    % ---------- Proposal sample v^(k)---------------------------------
    new_online_parameters = gen_proposal(offline_data);
%     new_online_parameters = ...
%                             [...
%                             10^( 9+3*rand() ), ...
%                             10^( 2+3*rand() ), ...
%                             2000+1000*rand(), ...
%                             5+45*rand(), ...
%                             10^( 3+4*rand() ) ...
%                             ];
    % ------------------------------------
    
    % Observation G(new_u)
    [new_observation, new_time_of_erupt] = MLMCMC_fun_observation_explicit_RK( offline_data, real_time_step, (MCMC_time_level-1), new_online_parameters );
    
    % Potential function \Phi^{J_l,l}(new_u;\delta)
    constants{1} = factor_of_real_time*real_time_step*( 2^(-(MCMC_time_level-1)) );
    [ new_Phi, new_time_to_erupt ] = MLMCMC_fun_potential( constants, delta, new_observation, new_time_of_erupt );
    % -----------------------------------------------------------------
    
    
    % ------------- Acceptance probability ------------------
    e = exp( old_Phi - new_Phi );
    alpha = min(1,e);
    % -----------------------------------------------------------------
    
    
    % -------------------- Decision -------------------------
    if rand <= alpha        %Accept the proposal (i.e. u^{k+1) = v^(k))
        old_online_parameters = new_online_parameters;      %<---- Note that we always update old_u to be the latest sample
        old_time_to_erupt = new_time_to_erupt;
        old_observation = new_observation;
        old_Phi = new_Phi;
        
        col_time_to_erupt_part(i) = old_time_to_erupt;
        col_final_tilt_part(i) = old_observation(end);
        
        % ********* Preparation for computing A_{2}^l(v^(k)) *****************
        %  +++++++ To compute \Phi^{J_l,l}(v^(k);\delta) ++++++
        % Observation G^{J_l,l}(old_u)
        [old_observation_p1, old_time_of_erupt_p1] = MLMCMC_fun_observation_explicit_RK( offline_data, real_time_step, MCMC_time_level, old_online_parameters );
        
        % Potential function \Phi^{J_l,l}(old_u;\delta)
        constants{1} = factor_of_real_time*real_time_step*( 2^(-MCMC_time_level) );
        [ old_Phi_p1, old_time_to_erupt_p1 ] = MLMCMC_fun_potential( constants, delta, old_observation_p1, old_time_of_erupt_p1 );
        % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % ****************************************************************
        
        % Note old_Phi_p1-old_Phi = \Phi^{J_l, l}(u;\delta) - \Phi^{J_{l-1}, l-1}(u;\delta)
        I_l = old_Phi_p1 - old_Phi;
        if  I_l > 0     %<------------------------- Different
            % A_{2}^{l}(v^(k)) = 1 - exp[\Phi^{J_{l-1}, l-1}(v^(k);\delta) - \Phi^{J_l, l}(v^(k);\delta)]
            A2_l = 1 - exp( -1*I_l );
            
            %  B^{l'}(v^(k))
            B = B_lprime_when_lprime_is_1( old_online_parameters, old_time_to_erupt, old_observation );
            
            % Collection of results
            cum_E_l1_A2_B_1nI = cum_E_l1_A2_B_1nI + A2_l*B;
            cum_E_l1_A2_1nI = cum_E_l1_A2_1nI + A2_l;
            cum_E_l1_1nA2_B_1nI = cum_E_l1_1nA2_B_1nI + (1 - A2_l)*B;
        else
            %  B^{l'}(v^(k))
            B = B_lprime_when_lprime_is_1( old_online_parameters, old_time_to_erupt, old_observation );
            cum_E_l1_B_I = cum_E_l1_B_I + B;
        end
    else
        % When the (i+1)th sample rejects the proposal sample, i.e. takes
        % the old sample hence need not compute anything new
        
        col_time_to_erupt_part(i) = old_time_to_erupt;
        col_final_tilt_part(i) = old_observation(end);
        
        if  I_l > 0     %<------------------------- Different
            cum_E_l1_A2_B_1nI = cum_E_l1_A2_B_1nI + A2_l*B;
            cum_E_l1_A2_1nI = cum_E_l1_A2_1nI + A2_l;
            cum_E_l1_1nA2_B_1nI = cum_E_l1_1nA2_B_1nI + (1 - A2_l)*B;
        else
            cum_E_l1_B_I = cum_E_l1_B_I + B;
        end
    end
    
    i = i + 1;
    
end

E_l1_A2_B_1nI = cum_E_l1_A2_B_1nI/M_ll;
E_l1_A2_1nI = cum_E_l1_A2_1nI/M_ll;
E_l1_1nA2_B_1nI = cum_E_l1_1nA2_B_1nI/M_ll;
E_l1_B_I = cum_E_l1_B_I/M_ll;
end

function add_on = c_ll_rho_lminus1_when_lprime_is_1(alpha, l, l_prime)
%here alpha = 0,1,2,3 or 4
% if alpha==4
%     add_on = l^2;
% else
%     add_on = (2-0.5*alpha)^(-2) * l^(alpha-2);
% end

% ----- one telescopic sum -----
% alpha | total factor error
%   0   |       O(L)
%   1   |       O(L^(1/2))
%   2   |       O(ln(L))
add_on = l^(alpha);
% ------------------------------

end

% ##################################################################################################################
%% ------ Plain_MCMC_rho_l_when_lprime_is_1 ---------------------
function [ E_l2_A1_B_I, E_l2_A1_I, E_l2_1nA1_B_I, E_l2_B_1nI, col_time_to_erupt_part, col_final_tilt_part ] = Plain_MCMC_rho_l_when_lprime_is_1( offline_data, real_time_step, MCMC_time_level, delta, L, rate, alpha, initial_old_online_parameters_wrt_l )
%Plain_MCMC_rho_l_when_lprime_is_1 Computes the estimated expectation of
%A_{1}^{l}(u)*B^{1}(u)*I^{l}(u), A_{1}^{l}(u)*I^{l}(u), 
%(1 - A_{1}^{l}(u))*B^{1}(u)*I^{l}(u), B^{1}(u)*( 1-I^{l}(u) )
%with respect to the posterior probability measure \rho^{J_l, l, \delta}
%   [L_nodes, L_weights] - nodes and weights of the Gauss Legendre quad
%   left, right - left and right endpoint of the domain
%   delta - data (Computed by 1 realization of u and \vartheta)
%   L - Levels of telescopic expansions (or the finest mesh-size 2^{-L})
%   l, l_prime - current level or term in the double summation of E_L^{MLMCMC}
%   rate - rate of convergence of Observation and Quantity of Interest (if
%          the convergence is O(h^r), then rate = r)
%
%   MCMC_time_level represents l from before


%       *Note l_prime == 1 here

% Collect results cummulatively
cum_E_l2_A1_B_I = 0;
cum_E_l2_A1_I = 0;
cum_E_l2_1nA1_B_I = 0;
cum_E_l2_B_1nI = 0;

l_prime = 1;

% Plain MCMC sampling size
M_ll = ceil( c_ll_rho_l_when_lprime_is_1(alpha, MCMC_time_level, l_prime)*2^( 2*rate*( L-MCMC_time_level ) ) );

col_time_to_erupt_part = zeros(M_ll,1);
col_final_tilt_part = zeros(M_ll,1);

% ------- Initializing constants -------
factor_of_real_time = 4;
MCMC_time_step = factor_of_real_time*real_time_step*(2^(-MCMC_time_level));

% Number of observation
no_of_obs = length(delta);

constants = {MCMC_time_step, real_time_step, no_of_obs, offline_data.constants{10} };
% --------------------------------------

% ------ Generate initial sample of the Markov Chain ------------------------------------------------------------
old_online_parameters = initial_old_online_parameters_wrt_l{1};
old_time_to_erupt = initial_old_online_parameters_wrt_l{2};
old_observation = initial_old_online_parameters_wrt_l{3};
old_Phi = initial_old_online_parameters_wrt_l{4};

    % ***** A_{1}^l(old_u) = 1 - exp[\Phi^{J_l, l}(old_u;\delta) - \Phi^{J_{l-1}, l-1}(old_u;\delta)] *****
    % To compute \Phi^{J_{l-1},l-1}(old_u;\delta)
    % Observation G^{J_{l-1},l-1}(old_u)
    [old_observation_n1, old_time_of_erupt_n1] = MLMCMC_fun_observation_explicit_RK( offline_data, real_time_step, (MCMC_time_level-1), old_online_parameters );

    % Potential function \Phi^{J_{l-1},l-1}(old_u;\delta)
    constants{1} = factor_of_real_time*real_time_step*( 2^( -(MCMC_time_level-1) ) );
    [ old_Phi_n1, old_time_to_erupt_n1 ] = MLMCMC_fun_potential( constants, delta, old_observation_n1, old_time_of_erupt_n1 );

    % Recall I^{l}(u) = { 1 if \Phi^{J_l, l}(u;\delta) - \Phi^{J_{l-1}, l-1}(u;\delta) <= 0
    %                   { 0 otherwise
    % This acts like a decision whether to use a sample a not, which will
    % be used later
    I_l = old_Phi - old_Phi_n1;
    
    % A_{1}^{l}(old_u) = 1 - exp[\Phi^{J_l, l}(old_u;\delta) - \Phi^{J_{l-1}, l-1}(old_u;\delta)]
    A1_l = 1 - exp( I_l );
    % ******************************************************************************************************

    % ########------ B^{1}(old_u) -----##########################
    B = B_lprime_when_lprime_is_1( old_online_parameters, old_time_to_erupt, old_observation );
    % ############################################################
% ---------------------------------------------------------------------------------------------------------------------

% -------- coefficient of v^(k) in pCN -----------------
% weight_old_u = sqrt(1-beta^2);
% -------------------------------------------------------
% no_online_para = size(initial_old_online_parameters_wrt_l);
i = 1;
while i <= M_ll
    % ---------- Proposal sample v^(k)---------------------------------
    new_online_parameters = gen_proposal(offline_data);
%     new_online_parameters = ...
%                             [...
%                             10^( 9+3*rand() ), ...
%                             10^( 2+3*rand() ), ...
%                             2000+1000*rand(), ...
%                             5+45*rand(), ...
%                             10^( 3+4*rand() ) ...
%                             ];
        
    % Observation G(new_u)
    [new_observation, new_time_of_erupt] = MLMCMC_fun_observation_explicit_RK( offline_data, real_time_step, MCMC_time_level, new_online_parameters );
    
    % Potential function \Phi^{J_l,l}(new_u;\delta)
    constants{1} = factor_of_real_time*real_time_step*( 2^(-MCMC_time_level) );
    [ new_Phi, new_time_to_erupt ] = MLMCMC_fun_potential( constants, delta, new_observation, new_time_of_erupt );
    % -----------------------------------------------------------------
    
    
    % ------------- Acceptance probability ------------------
    e = exp( old_Phi - new_Phi );
    alpha = min(1,e);
    % -----------------------------------------------------------------
    

    % -------------------- Decision -------------------------
    if rand <= alpha        %Accept the proposal (i.e. u^{k+1) = v^(k))
        old_online_parameters = new_online_parameters;      %<---- Note that we always update old_u to be the latest sample
        old_time_to_erupt = new_time_to_erupt;
        old_observation = new_observation;
        old_Phi = new_Phi;
        
        col_time_to_erupt_part(i) = old_time_to_erupt;
        col_final_tilt_part(i) = old_observation(end);
        
        % ********* Preparation for computing A^l(v^(k)) *****************
        %  +++++++ To compute \Phi^{J_{l-1},l-1}(v^(k);\delta) ++++++
        % Observation G^{J_{l-1},l-1}(old_u)
        [old_observation_n1, old_time_of_erupt_n1] = MLMCMC_fun_observation_explicit_RK( offline_data, real_time_step, (MCMC_time_level-1), old_online_parameters );

        % Potential function \Phi^{J_{l-1},l-1}(old_u;\delta)
        constants{1} = factor_of_real_time*real_time_step*( 2^(-(MCMC_time_level-1)) );
        [ old_Phi_n1, old_time_to_erupt_n1 ] = MLMCMC_fun_potential( constants, delta, old_observation_n1, old_time_of_erupt_n1 );
        % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % ****************************************************************
        
        I_l = old_Phi - old_Phi_n1;
        if  I_l <= 0
            % A_{1}^{l}(v^(k)) = 1 - exp[\Phi^{J_l, l}(v^(k);\delta) - \Phi^{J_{l-1}, l-1}(v^(k);\delta)]
            A1_l = 1 - exp( I_l );
            
            %  B^{1}(v^(k))
            B = B_lprime_when_lprime_is_1( old_online_parameters, old_time_to_erupt, old_observation );
            
            % Collection of results
            cum_E_l2_A1_B_I = cum_E_l2_A1_B_I + (A1_l * B);
            cum_E_l2_A1_I = cum_E_l2_A1_I + A1_l;
            cum_E_l2_1nA1_B_I = cum_E_l2_1nA1_B_I + ( (1 - A1_l) * B );
        else
            %  B^{1}(v^(k))
            B = B_lprime_when_lprime_is_1( old_online_parameters, old_time_to_erupt, old_observation );
            cum_E_l2_B_1nI = cum_E_l2_B_1nI + B;
        end
    else
        % When the (i+1)th sample rejects the proposal sample, i.e. takes
        % the old sample hence need not compute anything new
        col_time_to_erupt_part(i) = old_time_to_erupt;
        col_final_tilt_part(i) = old_observation(end);
        
        if  I_l <= 0
            cum_E_l2_A1_B_I = cum_E_l2_A1_B_I + (A1_l * B);
            cum_E_l2_A1_I = cum_E_l2_A1_I + A1_l;
            cum_E_l2_1nA1_B_I = cum_E_l2_1nA1_B_I + ( (1 - A1_l) * B );
        else
            cum_E_l2_B_1nI = cum_E_l2_B_1nI + B;
        end
    end
    
    i = i + 1;
    
end

E_l2_A1_B_I = cum_E_l2_A1_B_I/M_ll;
E_l2_A1_I = cum_E_l2_A1_I/M_ll;
E_l2_1nA1_B_I = cum_E_l2_1nA1_B_I/M_ll;
E_l2_B_1nI = cum_E_l2_B_1nI/M_ll;
end

function add_on = c_ll_rho_l_when_lprime_is_1(alpha, l, l_prime)
%here alpha = 0,1,2,3 or 4
% if alpha==4
%     add_on = l^2;
% else
%     add_on = (2-0.5*alpha)^(-2) * l^(alpha-2);
% end

% ----- one telescopic sum -----
% alpha | total factor error
%   0   |       O(L)
%   1   |       O(L^(1/2))
%   2   |       O(ln(L))
add_on = l^(alpha);
% ------------------------------
end

% ##################################################################################################################
%% ----- Plain_MCMC_B_wrt_rho_l_when_l_and_lprime_are_1 -----
function [ E_B, col_time_to_erupt_part, col_final_tilt_part ] = Plain_MCMC_B_wrt_rho_l_when_l_and_lprime_are_1( offline_data, real_time_step, MCMC_time_level, delta, L, rate, alpha, initial_old_online_parameters_wrt_l )
%Plain_MCMC_B_when_l_and_lprime_are_1 Computes the estimated expectation of B^{1}(u)
%with respect to the posterior probability measure \rho^{J_1, 1, \delta}
%   [L_nodes, L_weights] - nodes and weights of the Gauss Legendre quad
%   left, right - left and right endpoint of the domain
%   delta - data (Computed by 1 realization of u and \vartheta)
%   L - Levels of telescopic expansions (or the finest mesh-size 2^{-L})
%   l, l_prime - current level or term in the double summation of E_L^{MLMCMC}
%   rate - rate of convergence of Observation and Quantity of Interest(if
%          the convergence is O(h^r), then rate = r)
%
%   MCMC_time_level represents l from before

%       *The l in the input should always be 1, I left it as a input
%       variable for easy referencing later, like the sample size M_{ll'}
%           And also, l_prime = 1. 


% cum_E_B - collect results cummulatively
cum_E_B = 0;

% Plain MCMC sampling size
%M_ll = 2^( 2*rate*( L-( l+l_prime) ) );
M_ll = ceil( c_ll_rho_l_when_l_and_lprime_are_1(alpha,L) * 2^(2*rate*( L )) );  % <----- Proper sampling size

% ----- 12/02/2020 -----
col_time_to_erupt_part = zeros(M_ll,1);
col_final_tilt_part = zeros(M_ll,1);
% ---------------------

% ------- Initializing constants -------
factor_of_real_time = 4;
MCMC_time_step = factor_of_real_time*real_time_step*(2^(-MCMC_time_level));

% No of observation
no_of_obs = length(delta);

constants = {MCMC_time_step, real_time_step, no_of_obs, offline_data.constants{10} };
% --------------------------------------


% ------ Initial sample of the Markov Chain from Burn in------------
old_online_parameters = initial_old_online_parameters_wrt_l{1};
old_time_to_erupt = initial_old_online_parameters_wrt_l{2};
old_observation = initial_old_online_parameters_wrt_l{3};
old_Phi = initial_old_online_parameters_wrt_l{4};

% B^{l'}(old_u)
B = B_lprime_when_lprime_is_1( old_online_parameters, old_time_to_erupt, old_observation );
% --------------------------------------------------------------

% sqrt_1nBsq = sqrt(1-beta^2);
% no_online_para = size(initial_old_online_parameters_wrt_l);
i = 1;
while i <= M_ll   
    %  --- offline_data.online_parameters_property = { op_range, op_mean,  op_percent_range_to_move } ---
    new_online_parameters = gen_proposal(offline_data);
%     new_online_parameters = ...
%                             [...
%                             10^( 9+3*rand() ), ...
%                             10^( 2+3*rand() ), ...
%                             2000+1000*rand(), ...
%                             5+45*rand(), ...
%                             10^( 3+4*rand() ) ...
%                             ];
%     collect_online_para(i,:) = new_online_parameters;
    
    % Observation G(new_u)
    [new_observation, new_time_of_erupt] = MLMCMC_fun_observation_explicit_RK( offline_data, real_time_step, MCMC_time_level, new_online_parameters );

    % Potential function \Phi^{J_l,l}(new_u;\delta)
%     new_Phi = 0.5*(delta - new_observation)^2;
    
    [ new_Phi, new_time_to_erupt ] = MLMCMC_fun_potential( constants, delta, new_observation, new_time_of_erupt );
    % ------------------------------------------------------------------
    
    
    % ------------- Acceptance probability ------------------
    e = exp( old_Phi - new_Phi );
    alpha = min(1,e);
    %--------------------------------------------------------
    
    
    % ------------ Decision ----------------
    if rand <= alpha    %Accept the proposal (i.e. u^{k+1) = v^(k))
        old_online_parameters = new_online_parameters;
        old_time_to_erupt = new_time_to_erupt;
        old_observation = new_observation;
        old_Phi = new_Phi;
        
        col_time_to_erupt_part(i) = old_time_to_erupt;
        col_final_tilt_part(i) = old_observation(end);
        
        B = B_lprime_when_lprime_is_1( old_online_parameters, old_time_to_erupt, old_observation );
        cum_E_B = cum_E_B + B;
    else    %Reject the proposal(use back the old sample; i.e. u^{k+1) = u^(k))
        col_time_to_erupt_part(i) = old_time_to_erupt;
        col_final_tilt_part(i) = old_observation(end);
        
        cum_E_B = cum_E_B + B;
    end
    % --------------------------------------
    i = i + 1;
end

E_B = cum_E_B/M_ll;
end

function add_on = c_ll_rho_l_when_l_and_lprime_are_1(alpha, L)
%here we only use alpha = 1,2,3 or 4
% if alpha == 4
%     add_on = 1/(log(L))^2;
% else
%     add_on = L^(alpha-4);
% end

% ----- one telescopic sum -----
add_on = (alpha==0)*L^(-2) + (alpha==1)*L^(-1) + (alpha==2)*( log(L)^(-2) );
% ------------------------------
end
% ##################################################################################################################
%% ------ Burn_in_lognormal ------------
function [ initial_parameters ] = Burn_in_lognormal( offline_data, online_parameters, real_time_step, MCMC_time_level, delta, burn )
%UNTITLED Summary of this function goes here
%   offline_data.constants = {P_atm, g, nu, dist, Rgas, Temperature, InitialMagmaLevel, Ls0, order_of_ODE, no_of_RKGL_nodes, obs_cov_mat}
%   offline_data.quadrature = {fun_list, dfun_list, RK_A, RK_b, RK_c}
%   offline_data.online_parameters_property = { op_range, op_mean, op_percent_range_to_move }
%   online_parameters = {G, mu, rho, mass0, rc}
%   real_time_step must be a multiple or factor of MCMC_time_step


% ------- Initializing constants -------
factor_of_real_time = 4;
MCMC_time_step = factor_of_real_time*real_time_step*(2^(-MCMC_time_level));
% No of observation
no_of_obs = length(delta);
constants = {MCMC_time_step, real_time_step, no_of_obs, offline_data.constants{10} };
% --------------------------------------

% obs_cov_mat = offline_data.constants{10};



% ------ Generate initial sample of the Markov Chain ------------
old_online_parameters = online_parameters;


% Observation G(old_u)
[old_observation, old_time_of_erupt] = MLMCMC_fun_observation_explicit_RK( offline_data, real_time_step, MCMC_time_level, old_online_parameters );

% --- Potential function \Phi^{J_l,l}(old_u;\delta) ---
% old_Phi = 0.5*(delta - old_observation)^2;
    %   constants - {MCMC_time_step, real_time_step, no_of_obs, obs_cov_mat }
    %   obs_cov_mat - offline_data.constants{10};
    [ old_Phi, old_time_to_erupt ] = MLMCMC_fun_potential( constants, delta, old_observation, old_time_of_erupt );

% sqrt_1nBsq = sqrt(1-beta^2);

% ------------------- Burn in ------------------------------
% collect_online_para = zeros(burn,length(online_parameters));
% collect_MCMC_para = zeros(burn,length(online_parameters));
% no_online_para = size(online_parameters);

for i = 1:burn  
    % ------ pCN v^(k) -------------------
%     new_u = sqrt_1nBsq*old_u + beta*randn;  %replace randn by inv(prior_cov_mat).^(1/2) * randn
    % ------------------------------------
    
    
    %  --- offline_data.online_parameters_property = { op_range, op_mean,  op_percent_range_to_move } ---
    new_online_parameters = gen_proposal(offline_data);
%     new_online_parameters = ...
%                             [...
%                             10^( 9+3*rand() ), ...
%                             10^( 2+3*rand() ), ...
%                             2000+1000*rand(), ...
%                             5+45*rand(), ...
%                             10^( 3+4*rand() ) ...
%                             ];
%     collect_online_para(i,:) = new_online_parameters;

    
    % Observation G(new_u)
    [new_observation, new_time_of_erupt] = MLMCMC_fun_observation_explicit_RK( offline_data, real_time_step, MCMC_time_level, new_online_parameters );

    % Potential function \Phi^{J_l,l}(new_u;\delta)
%     new_Phi = 0.5*(delta - new_observation)^2;
    
    [ new_Phi, new_time_to_erupt ] = MLMCMC_fun_potential( constants, delta, new_observation, new_time_of_erupt );
    % ------------------------------------------------------------------


    % ------------- Acceptance probability ------------------
    e = exp( old_Phi - new_Phi );
    alpha = min(1,e);
    % -------------------------------------------------------


    % -------------- Decision --------------
    if rand <= alpha    %Accept the proposal (i.e. u^{k+1) = v^(k))
        old_online_parameters = new_online_parameters;  %<---- Note that we always update old_online_parameters to be the latest sample
        old_time_to_erupt = new_time_to_erupt;
        old_observation = new_observation;
        old_Phi = new_Phi;
    end
%     collect_MCMC_para(i,:)=old_online_parameters;
    % ------------------------------------
    
end
initial_parameters = {old_online_parameters, old_time_to_erupt, old_observation, old_Phi};
end
% ##################################################################################################################

% ##################################################################################################################
% ----- function to generate proposal -----------
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

% ##################################################################################################################
function [approx_index_of_roots] = root_finder(delta, full_observation, diff_mean40_mean60)
    global trigger_condition
    global standard_dev
    if diff_mean40_mean60 < trigger_condition
%         full_observation_minus_delta_end = full_observation - mean( delta( (end-floor(0.1*length(delta))) :end) );
        full_observation_minus_delta_end = full_observation - mean( delta( (end-10) :end) );
    else
        full_observation_minus_delta_end = full_observation - delta(end);
    end
    
    len_full_obs = length(full_observation);
    [ min_value, min_point ] = min( full_observation );
    no_of_roots = (full_observation_minus_delta_end(1)*full_observation_minus_delta_end(min_point)<0)...
                    +(full_observation_minus_delta_end(min_point)*full_observation_minus_delta_end(end)<0);
    
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
        [min_diff, approx_index_of_roots] = min(abs(full_observation-delta(end)));
    else
        error('There are more than 2 roots!')
    end
end

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

% function [ min_point ] = min_point_finder( full_observation )
% %UNTITLED Summary of this function goes here
% %   Assumption: tilt is either increasing or 
% %               slightly decreasing at first and then increasing
% 
%     if full_observation(end) - full_observation(1)<0
%         error('Last observation is lesser than first observation!')
%     end
% 
%     grad_obs = full_observation(2:end) - full_observation(1:(end-1));
%     if grad_obs(1) > 0   %tilt is increasing
%         min_point = 1;
%     else
%         min_point = grad_finder(grad_obs, 1, length(grad_obs) );
%     end
% 
% end
% ##################################################################################################################
