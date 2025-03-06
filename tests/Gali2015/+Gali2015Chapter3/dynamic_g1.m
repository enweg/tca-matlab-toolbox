function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = Gali2015Chapter3.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(25, 36);
g1(1,7)=1;
g1(1,32)=(-params(2));
g1(1,8)=(-T(4));
g1(2,32)=T(3);
g1(2,8)=1;
g1(2,33)=(-1);
g1(2,12)=T(3);
g1(2,14)=(-T(3));
g1(3,7)=(-params(8));
g1(3,11)=(-params(9));
g1(3,14)=1;
g1(3,19)=(-1);
g1(4,12)=1;
g1(4,20)=(-((1-params(3))*T(1)*(-params(6))));
g1(4,25)=(-(1-params(5)));
g1(5,32)=1;
g1(5,13)=1;
g1(5,14)=(-1);
g1(6,9)=1;
g1(6,20)=(-T(1));
g1(7,8)=1;
g1(7,9)=1;
g1(7,10)=(-1);
g1(8,3)=(-params(4));
g1(8,19)=1;
g1(8,35)=(-1);
g1(9,4)=(-params(3));
g1(9,20)=1;
g1(9,34)=(-1);
g1(10,10)=1;
g1(10,15)=(-(1-params(1)));
g1(10,20)=(-1);
g1(11,5)=(-params(5));
g1(11,25)=1;
g1(11,36)=1;
g1(12,7)=(-4);
g1(12,1)=4;
g1(12,10)=(-4);
g1(12,2)=(-(4*params(10)));
g1(12,14)=(-(4*(-params(10))));
g1(12,17)=1;
g1(13,10)=(-1);
g1(13,14)=params(10);
g1(13,16)=1;
g1(14,14)=(-4);
g1(14,22)=1;
g1(15,13)=(-4);
g1(15,21)=1;
g1(16,12)=(-4);
g1(16,23)=1;
g1(17,7)=(-4);
g1(17,24)=1;
g1(18,10)=(-1);
g1(18,11)=1;
g1(19,7)=1;
g1(19,6)=1;
g1(19,26)=(-1);
g1(20,10)=1;
g1(20,28)=(-1);
g1(21,15)=(-params(7));
g1(21,26)=(-1);
g1(21,27)=1;
g1(21,28)=(-params(6));
g1(22,26)=1;
g1(22,27)=(-1);
g1(22,29)=1;
g1(23,16)=(-1);
g1(23,18)=1;
g1(23,26)=(-1);
g1(24,10)=T(2);
g1(24,20)=(-((1+params(7))/(1-params(1))));
g1(24,30)=1;
g1(25,8)=T(2);
g1(25,31)=1;

end
