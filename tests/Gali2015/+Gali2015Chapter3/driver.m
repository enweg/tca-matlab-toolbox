%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'Gali2015Chapter3';
M_.dynare_version = '5.5-arm64';
oo_.dynare_version = '5.5-arm64';
options_.dynare_version = '5.5-arm64';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(3,1);
M_.exo_names_tex = cell(3,1);
M_.exo_names_long = cell(3,1);
M_.exo_names(1) = {'eps_a'};
M_.exo_names_tex(1) = {'{\varepsilon_a}'};
M_.exo_names_long(1) = {'technology shock'};
M_.exo_names(2) = {'eps_nu'};
M_.exo_names_tex(2) = {'{\varepsilon_\nu}'};
M_.exo_names_long(2) = {'monetary policy shock'};
M_.exo_names(3) = {'eps_z'};
M_.exo_names_tex(3) = {'{\varepsilon_z}'};
M_.exo_names_long(3) = {'preference shock innovation'};
M_.endo_names = cell(25,1);
M_.endo_names_tex = cell(25,1);
M_.endo_names_long = cell(25,1);
M_.endo_names(1) = {'pi'};
M_.endo_names_tex(1) = {'{\pi}'};
M_.endo_names_long(1) = {'inflation'};
M_.endo_names(2) = {'y_gap'};
M_.endo_names_tex(2) = {'{\tilde y}'};
M_.endo_names_long(2) = {'output gap'};
M_.endo_names(3) = {'y_nat'};
M_.endo_names_tex(3) = {'{y^{nat}}'};
M_.endo_names_long(3) = {'natural output'};
M_.endo_names(4) = {'y'};
M_.endo_names_tex(4) = {'{y}'};
M_.endo_names_long(4) = {'output'};
M_.endo_names(5) = {'yhat'};
M_.endo_names_tex(5) = {'{\hat y}'};
M_.endo_names_long(5) = {'output deviation from steady state'};
M_.endo_names(6) = {'r_nat'};
M_.endo_names_tex(6) = {'{r^{nat}}'};
M_.endo_names_long(6) = {'natural interest rate'};
M_.endo_names(7) = {'r_real'};
M_.endo_names_tex(7) = {'{r^r}'};
M_.endo_names_long(7) = {'real interest rate'};
M_.endo_names(8) = {'i'};
M_.endo_names_tex(8) = {'{i}'};
M_.endo_names_long(8) = {'nominal interest rate'};
M_.endo_names(9) = {'n'};
M_.endo_names_tex(9) = {'{n}'};
M_.endo_names_long(9) = {'hours worked'};
M_.endo_names(10) = {'m_real'};
M_.endo_names_tex(10) = {'{m-p}'};
M_.endo_names_long(10) = {'real money stock'};
M_.endo_names(11) = {'m_growth_ann'};
M_.endo_names_tex(11) = {'{\Delta m}'};
M_.endo_names_long(11) = {'money growth annualized'};
M_.endo_names(12) = {'m_nominal'};
M_.endo_names_tex(12) = {'{m}'};
M_.endo_names_long(12) = {'nominal money stock'};
M_.endo_names(13) = {'nu'};
M_.endo_names_tex(13) = {'{\nu}'};
M_.endo_names_long(13) = {'AR(1) monetary policy shock process'};
M_.endo_names(14) = {'a'};
M_.endo_names_tex(14) = {'{a}'};
M_.endo_names_long(14) = {'AR(1) technology shock process'};
M_.endo_names(15) = {'r_real_ann'};
M_.endo_names_tex(15) = {'{r^{r,ann}}'};
M_.endo_names_long(15) = {'annualized real interest rate'};
M_.endo_names(16) = {'i_ann'};
M_.endo_names_tex(16) = {'{i^{ann}}'};
M_.endo_names_long(16) = {'annualized nominal interest rate'};
M_.endo_names(17) = {'r_nat_ann'};
M_.endo_names_tex(17) = {'{r^{nat,ann}}'};
M_.endo_names_long(17) = {'annualized natural interest rate'};
M_.endo_names(18) = {'pi_ann'};
M_.endo_names_tex(18) = {'{\pi^{ann}}'};
M_.endo_names_long(18) = {'annualized inflation rate'};
M_.endo_names(19) = {'z'};
M_.endo_names_tex(19) = {'{z}'};
M_.endo_names_long(19) = {'AR(1) preference shock process'};
M_.endo_names(20) = {'p'};
M_.endo_names_tex(20) = {'{p}'};
M_.endo_names_long(20) = {'price level'};
M_.endo_names(21) = {'w'};
M_.endo_names_tex(21) = {'{w}'};
M_.endo_names_long(21) = {'nominal wage'};
M_.endo_names(22) = {'c'};
M_.endo_names_tex(22) = {'{c}'};
M_.endo_names_long(22) = {'consumption'};
M_.endo_names(23) = {'w_real'};
M_.endo_names_tex(23) = {'{\frac{w}{p}}'};
M_.endo_names_long(23) = {'real wage'};
M_.endo_names(24) = {'mu'};
M_.endo_names_tex(24) = {'{\mu}'};
M_.endo_names_long(24) = {'markup'};
M_.endo_names(25) = {'mu_hat'};
M_.endo_names_tex(25) = {'{\hat \mu}'};
M_.endo_names_long(25) = {'markup gap'};
M_.endo_partitions = struct();
M_.param_names = cell(12,1);
M_.param_names_tex = cell(12,1);
M_.param_names_long = cell(12,1);
M_.param_names(1) = {'alppha'};
M_.param_names_tex(1) = {'{\alpha}'};
M_.param_names_long(1) = {'capital share'};
M_.param_names(2) = {'betta'};
M_.param_names_tex(2) = {'{\beta}'};
M_.param_names_long(2) = {'discount factor'};
M_.param_names(3) = {'rho_a'};
M_.param_names_tex(3) = {'{\rho_a}'};
M_.param_names_long(3) = {'autocorrelation technology shock'};
M_.param_names(4) = {'rho_nu'};
M_.param_names_tex(4) = {'{\rho_{\nu}}'};
M_.param_names_long(4) = {'autocorrelation monetary policy shock'};
M_.param_names(5) = {'rho_z'};
M_.param_names_tex(5) = {'{\rho_{z}}'};
M_.param_names_long(5) = {'autocorrelation preference shock'};
M_.param_names(6) = {'siggma'};
M_.param_names_tex(6) = {'{\sigma}'};
M_.param_names_long(6) = {'inverse EIS'};
M_.param_names(7) = {'varphi'};
M_.param_names_tex(7) = {'{\varphi}'};
M_.param_names_long(7) = {'inverse Frisch elasticity'};
M_.param_names(8) = {'phi_pi'};
M_.param_names_tex(8) = {'{\phi_{\pi}}'};
M_.param_names_long(8) = {'inflation feedback Taylor Rule'};
M_.param_names(9) = {'phi_y'};
M_.param_names_tex(9) = {'{\phi_{y}}'};
M_.param_names_long(9) = {'output feedback Taylor Rule'};
M_.param_names(10) = {'eta'};
M_.param_names_tex(10) = {'{\eta}'};
M_.param_names_long(10) = {'semi-elasticity of money demand'};
M_.param_names(11) = {'epsilon'};
M_.param_names_tex(11) = {'{\epsilon}'};
M_.param_names_long(11) = {'demand elasticity'};
M_.param_names(12) = {'theta'};
M_.param_names_tex(12) = {'{\theta}'};
M_.param_names_long(12) = {'Calvo parameter'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 25;
M_.param_nbr = 12;
M_.orig_endo_nbr = 25;
M_.aux_vars = [];
options_.varobs = cell(3, 1);
options_.varobs(1)  = {'pi'};
options_.varobs(2)  = {'i'};
options_.varobs(3)  = {'y_gap'};
options_.varobs_id = [ 1 8 2  ];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(3, 3);
M_.Correlation_matrix = eye(3, 3);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = true;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.nonzero_hessian_eqs = [];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 25;
M_.eq_nbr = 25;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 7 32;
 0 8 33;
 0 9 0;
 1 10 0;
 0 11 0;
 0 12 0;
 0 13 0;
 2 14 0;
 0 15 0;
 0 16 0;
 0 17 0;
 0 18 0;
 3 19 0;
 4 20 0;
 0 21 0;
 0 22 0;
 0 23 0;
 0 24 0;
 5 25 0;
 6 26 0;
 0 27 0;
 0 28 0;
 0 29 0;
 0 30 0;
 0 31 0;]';
M_.nstatic = 17;
M_.nfwrd   = 2;
M_.npred   = 6;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 6;
M_.ndynamic   = 8;
M_.dynamic_tmp_nbr = [4; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'New Keynesian Phillips Curve eq. (22)' ;
  2 , 'name' , 'Dynamic IS Curve eq. (23)' ;
  3 , 'name' , 'Interest Rate Rule eq. (26)' ;
  4 , 'name' , 'Definition natural rate of interest eq. (24)' ;
  5 , 'name' , 'Definition real interest rate' ;
  6 , 'name' , 'Definition natural output, eq. (20)' ;
  7 , 'name' , 'Definition output gap' ;
  8 , 'name' , 'Monetary policy shock' ;
  9 , 'name' , 'TFP shock' ;
  10 , 'name' , 'Production function (eq. 14)' ;
  11 , 'name' , 'Preference shock, p. 54' ;
  12 , 'name' , 'Money growth (derived from eq. (4))' ;
  13 , 'name' , 'Real money demand (eq. 4)' ;
  14 , 'name' , 'Annualized nominal interest rate' ;
  15 , 'name' , 'Annualized real interest rate' ;
  16 , 'name' , 'Annualized natural interest rate' ;
  17 , 'name' , 'Annualized inflation' ;
  18 , 'name' , 'Output deviation from steady state' ;
  19 , 'name' , 'Definition price level' ;
  20 , 'name' , 'resource constraint, eq. (12)' ;
  21 , 'name' , 'FOC labor, eq. (2)' ;
  22 , 'name' , 'definition real wage' ;
  23 , 'name' , 'definition nominal money stock' ;
  24 , 'name' , 'average price markup, eq. (18)' ;
  25 , 'name' , 'average price markup, eq. (20)' ;
};
M_.mapping.pi.eqidx = [1 2 3 5 12 17 19 ];
M_.mapping.y_gap.eqidx = [1 2 7 25 ];
M_.mapping.y_nat.eqidx = [6 7 ];
M_.mapping.y.eqidx = [7 10 12 13 18 20 24 ];
M_.mapping.yhat.eqidx = [3 18 ];
M_.mapping.r_nat.eqidx = [2 4 16 ];
M_.mapping.r_real.eqidx = [5 15 ];
M_.mapping.i.eqidx = [2 3 5 12 13 14 ];
M_.mapping.n.eqidx = [10 21 ];
M_.mapping.m_real.eqidx = [13 23 ];
M_.mapping.m_growth_ann.eqidx = [12 ];
M_.mapping.m_nominal.eqidx = [23 ];
M_.mapping.nu.eqidx = [3 8 ];
M_.mapping.a.eqidx = [4 6 9 10 24 ];
M_.mapping.r_real_ann.eqidx = [15 ];
M_.mapping.i_ann.eqidx = [14 ];
M_.mapping.r_nat_ann.eqidx = [16 ];
M_.mapping.pi_ann.eqidx = [17 ];
M_.mapping.z.eqidx = [4 11 ];
M_.mapping.p.eqidx = [19 21 22 23 ];
M_.mapping.w.eqidx = [21 22 ];
M_.mapping.c.eqidx = [20 21 ];
M_.mapping.w_real.eqidx = [22 ];
M_.mapping.mu.eqidx = [24 ];
M_.mapping.mu_hat.eqidx = [25 ];
M_.mapping.eps_a.eqidx = [9 ];
M_.mapping.eps_nu.eqidx = [8 ];
M_.mapping.eps_z.eqidx = [11 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [4 8 13 14 19 20 ];
M_.exo_names_orig_ord = [1:3];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(25, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(12, 1);
M_.endo_trends = struct('deflator', cell(25, 1), 'log_deflator', cell(25, 1), 'growth_factor', cell(25, 1), 'log_growth_factor', cell(25, 1));
M_.NNZDerivatives = [74; 0; -1; ];
M_.static_tmp_nbr = [4; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(6) = 1;
siggma = M_.params(6);
M_.params(7) = 5;
varphi = M_.params(7);
M_.params(8) = 1.5;
phi_pi = M_.params(8);
M_.params(9) = 0.125;
phi_y = M_.params(9);
M_.params(12) = 0.75;
theta = M_.params(12);
M_.params(4) = 0.5;
rho_nu = M_.params(4);
M_.params(5) = 0.5;
rho_z = M_.params(5);
M_.params(3) = 0.9;
rho_a = M_.params(3);
M_.params(2) = 0.99;
betta = M_.params(2);
M_.params(10) = 3.77;
eta = M_.params(10);
M_.params(1) = 0.25;
alppha = M_.params(1);
M_.params(11) = 9;
epsilon = M_.params(11);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
M_.Sigma_e(2, 2) = 0.0625;
M_.Sigma_e(3, 3) = 0.25;
resid;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 20;
options_.order = 1;
var_list_ = {'pi';'i';'y_gap'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'Gali2015Chapter3_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Gali2015Chapter3_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Gali2015Chapter3_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Gali2015Chapter3_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Gali2015Chapter3_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Gali2015Chapter3_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Gali2015Chapter3_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
