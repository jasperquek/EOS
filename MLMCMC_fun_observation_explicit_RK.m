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