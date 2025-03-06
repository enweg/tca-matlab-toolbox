function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
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
%   g1
%

if T_flag
    T = Gali2015Chapter3.static_g1_tt(T, y, x, params);
end
g1 = zeros(25, 25);
g1(1,1)=1-params(2);
g1(1,2)=(-T(3));
g1(2,1)=T(4);
g1(2,6)=T(4);
g1(2,8)=(-T(4));
g1(3,1)=(-params(8));
g1(3,5)=(-params(9));
g1(3,8)=1;
g1(3,13)=(-1);
g1(4,6)=1;
g1(4,14)=(-((1-params(3))*T(1)*(-params(6))));
g1(4,19)=(-(1-params(5)));
g1(5,1)=1;
g1(5,7)=1;
g1(5,8)=(-1);
g1(6,3)=1;
g1(6,14)=(-T(1));
g1(7,2)=1;
g1(7,3)=1;
g1(7,4)=(-1);
g1(8,13)=1-params(4);
g1(9,14)=1-params(3);
g1(10,4)=1;
g1(10,9)=(-(1-params(1)));
g1(10,14)=(-1);
g1(11,19)=1-params(5);
g1(12,1)=(-4);
g1(12,11)=1;
g1(13,4)=(-1);
g1(13,8)=params(10);
g1(13,10)=1;
g1(14,8)=(-4);
g1(14,16)=1;
g1(15,7)=(-4);
g1(15,15)=1;
g1(16,6)=(-4);
g1(16,17)=1;
g1(17,1)=(-4);
g1(17,18)=1;
g1(18,5)=1;
g1(19,1)=1;
g1(20,4)=1;
g1(20,22)=(-1);
g1(21,9)=(-params(7));
g1(21,20)=(-1);
g1(21,21)=1;
g1(21,22)=(-params(6));
g1(22,20)=1;
g1(22,21)=(-1);
g1(22,23)=1;
g1(23,10)=(-1);
g1(23,12)=1;
g1(23,20)=(-1);
g1(24,4)=T(2);
g1(24,14)=(-((1+params(7))/(1-params(1))));
g1(24,24)=1;
g1(25,2)=T(2);
g1(25,25)=1;
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
