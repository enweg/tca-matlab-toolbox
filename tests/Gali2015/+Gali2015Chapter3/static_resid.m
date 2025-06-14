function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = Gali2015Chapter3.static_resid_tt(T, y, x, params);
end
residual = zeros(25, 1);
lhs = y(1);
rhs = params(2)*y(1)+T(3)*y(2);
residual(1) = lhs - rhs;
lhs = y(2);
rhs = y(2)+T(4)*(y(8)-y(1)-y(6));
residual(2) = lhs - rhs;
lhs = y(8);
rhs = y(1)*params(8)+params(9)*y(5)+y(13);
residual(3) = lhs - rhs;
lhs = y(6);
rhs = (1-params(5))*y(19)+y(14)*(1-params(3))*T(1)*(-params(6));
residual(4) = lhs - rhs;
lhs = y(7);
rhs = y(8)-y(1);
residual(5) = lhs - rhs;
lhs = y(3);
rhs = T(1)*y(14);
residual(6) = lhs - rhs;
lhs = y(2);
rhs = y(4)-y(3);
residual(7) = lhs - rhs;
lhs = y(13);
rhs = y(13)*params(4)+x(2);
residual(8) = lhs - rhs;
lhs = y(14);
rhs = y(14)*params(3)+x(1);
residual(9) = lhs - rhs;
lhs = y(4);
rhs = y(14)+(1-params(1))*y(9);
residual(10) = lhs - rhs;
lhs = y(19);
rhs = params(5)*y(19)-x(3);
residual(11) = lhs - rhs;
lhs = y(11);
rhs = y(1)*4;
residual(12) = lhs - rhs;
lhs = y(10);
rhs = y(4)-y(8)*params(10);
residual(13) = lhs - rhs;
lhs = y(16);
rhs = y(8)*4;
residual(14) = lhs - rhs;
lhs = y(15);
rhs = y(7)*4;
residual(15) = lhs - rhs;
lhs = y(17);
rhs = y(6)*4;
residual(16) = lhs - rhs;
lhs = y(18);
rhs = y(1)*4;
residual(17) = lhs - rhs;
lhs = y(5);
rhs = y(4)-(y(4));
residual(18) = lhs - rhs;
residual(19) = y(1);
lhs = y(4);
rhs = y(22);
residual(20) = lhs - rhs;
lhs = y(21)-y(20);
rhs = params(6)*y(22)+params(7)*y(9);
residual(21) = lhs - rhs;
lhs = y(23);
rhs = y(21)-y(20);
residual(22) = lhs - rhs;
lhs = y(12);
rhs = y(10)+y(20);
residual(23) = lhs - rhs;
lhs = y(24);
rhs = y(4)*(-T(2))+y(14)*(1+params(7))/(1-params(1));
residual(24) = lhs - rhs;
lhs = y(25);
rhs = y(2)*(-T(2));
residual(25) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
