function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = Gali2015Chapter3.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(25, 1);
lhs = y(7);
rhs = params(2)*y(32)+y(8)*T(4);
residual(1) = lhs - rhs;
lhs = y(8);
rhs = T(3)*(y(14)-y(32)-y(12))+y(33);
residual(2) = lhs - rhs;
lhs = y(14);
rhs = y(7)*params(8)+params(9)*y(11)+y(19);
residual(3) = lhs - rhs;
lhs = y(12);
rhs = (1-params(5))*y(25)+y(20)*(1-params(3))*T(1)*(-params(6));
residual(4) = lhs - rhs;
lhs = y(13);
rhs = y(14)-y(32);
residual(5) = lhs - rhs;
lhs = y(9);
rhs = T(1)*y(20);
residual(6) = lhs - rhs;
lhs = y(8);
rhs = y(10)-y(9);
residual(7) = lhs - rhs;
lhs = y(19);
rhs = params(4)*y(3)+x(it_, 2);
residual(8) = lhs - rhs;
lhs = y(20);
rhs = params(3)*y(4)+x(it_, 1);
residual(9) = lhs - rhs;
lhs = y(10);
rhs = y(20)+(1-params(1))*y(15);
residual(10) = lhs - rhs;
lhs = y(25);
rhs = params(5)*y(5)-x(it_, 3);
residual(11) = lhs - rhs;
lhs = y(17);
rhs = 4*(y(7)+y(10)-y(1)-params(10)*(y(14)-y(2)));
residual(12) = lhs - rhs;
lhs = y(16);
rhs = y(10)-y(14)*params(10);
residual(13) = lhs - rhs;
lhs = y(22);
rhs = y(14)*4;
residual(14) = lhs - rhs;
lhs = y(21);
rhs = y(13)*4;
residual(15) = lhs - rhs;
lhs = y(23);
rhs = y(12)*4;
residual(16) = lhs - rhs;
lhs = y(24);
rhs = y(7)*4;
residual(17) = lhs - rhs;
lhs = y(11);
rhs = y(10)-(steady_state(4));
residual(18) = lhs - rhs;
lhs = y(7);
rhs = y(26)-y(6);
residual(19) = lhs - rhs;
lhs = y(10);
rhs = y(28);
residual(20) = lhs - rhs;
lhs = y(27)-y(26);
rhs = params(6)*y(28)+params(7)*y(15);
residual(21) = lhs - rhs;
lhs = y(29);
rhs = y(27)-y(26);
residual(22) = lhs - rhs;
lhs = y(18);
rhs = y(16)+y(26);
residual(23) = lhs - rhs;
lhs = y(30);
rhs = y(10)*(-T(2))+y(20)*(1+params(7))/(1-params(1));
residual(24) = lhs - rhs;
lhs = y(31);
rhs = y(8)*(-T(2));
residual(25) = lhs - rhs;

end
