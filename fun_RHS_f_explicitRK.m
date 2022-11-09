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