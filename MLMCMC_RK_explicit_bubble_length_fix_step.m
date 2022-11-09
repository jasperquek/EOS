    function [ L_collect, h_collect, s_collect, t_collect ] = MLMCMC_RK_explicit_bubble_length_fix_step( order_of_ODE, fun_list, t_0, time_step, iteration_limit, L_0, d1, d2, d3, d4, d_h1, Vb, A1, h_min_height_of_magma_before_erupt )
    %UNTITLED2 Summary of this function goes here
    %   fun_list -  cell array containing function handles
    %   d1, d2, d3, d4, d_h1, u_s, A - constants
    %   t_0 - initial starting time
    %   time_step - time step
    %   L_0 - initial condition (L(t_0); L'(t_0))

    % [ fun_list ] = fun_RHS_f( order_of_ODE );
    % [ dfun_list ] = dfun_RHS_f_dLj( order_of_ODE );
    % [ RK_A, RK_b, RK_c ] = Runge_Kutta_butchertableau( no_of_RKGL_nodes );

    % no_iteration = (final_time - t_0)/time_step;
    
    fun1 = fun_list{1};
    fun2 = fun_list{2};

    L_collect = zeros(iteration_limit+1, order_of_ODE);
    h_collect = zeros(iteration_limit+1, 1);
    s_collect = zeros(iteration_limit+1, 1);
    t_collect = zeros(iteration_limit+1, 1);

    h = d_h1 - A1*L_0(1) - Vb*t_0;

    % ---- storing initial values of L, dL/dt, h, s-----
    L_collect(1,:) = L_0';
    h_collect(1) = h;
    s_collect(1) = 0;
    t_collect(1) = t_0;
    % ---------------------------------------------------

    L_n = L_0;
    for i = 1:iteration_limit
        t_n = t_0 + (i-1)*time_step;
        
        k1 = time_step*[fun1(d1, d2, d3, d4, d_h1, Vb, A1, t_n, L_n); fun2(d1, d2, d3, d4, d_h1, Vb, A1, t_n, L_n)];
        k2 = time_step*[fun1(d1, d2, d3, d4, d_h1, Vb, A1, (t_n + (time_step/2)), (L_n + (k1/2))); fun2(d1, d2, d3, d4, d_h1, Vb, A1, (t_n + (time_step/2)), (L_n + (k1/2)))];
        k3 = time_step*[fun1(d1, d2, d3, d4, d_h1, Vb, A1, (t_n + (time_step/2)), (L_n + (k2/2))); fun2(d1, d2, d3, d4, d_h1, Vb, A1, (t_n + (time_step/2)), (L_n + (k2/2)))];
        k4 = time_step*[fun1(d1, d2, d3, d4, d_h1, Vb, A1, (t_n + time_step), (L_n + k3)); fun2(d1, d2, d3, d4, d_h1, Vb, A1, (t_n + time_step), (L_n + k3))];
        L_np1 = L_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
        t_np1 = t_0 + (i)*time_step;
        
        
        L_collect(i+1,:) = L_np1';
        L_n = L_np1;

        h = d_h1 - A1*L_n(1) - Vb*t_n;
        h_collect(i+1) = h;
        s_collect(i+1) = Vb*t_n;
        t_collect(i+1) = t_np1;

        if h < h_min_height_of_magma_before_erupt
            break
        end


    end
    % L_final = L_np1;
    L_collect = L_collect(1:(i+1),:);
    h_collect = h_collect(1:(i+1),:);
    s_collect = s_collect(1:(i+1),:);
    t_collect = t_collect(1:(i+1),:);
    end