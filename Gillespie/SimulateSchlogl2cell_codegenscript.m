% SIMULATESCHLOGL2CELL_SCRIPT   Generate MEX-function SimulateSchlogl2cell_mex
%  from SimulateSchlogl2cell.
% 
% Script generated from project 'SimulateSchlogl2cell.prj' on 01-Dec-2017.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

%% Define argument types for entry-point 'SimulateSchlogl2cell'.
ARGS = cell(1,1);
ARGS{1} = cell(6,1);
ARGS{1}{1} = coder.typeof(0,[1 2]);
ARGS{1}{2} = coder.typeof(0,[1 2]);
ARGS_1_3 = struct;
ARGS_1_3.s = coder.typeof(0);
ARGS_1_3.K2 = coder.typeof(0);
ARGS_1_3.K = coder.typeof(0);
ARGS_1_3.a = coder.typeof(0);
ARGS_1_3.N = coder.typeof(0);
ARGS_1_3.g = coder.typeof(0);
ARGS{1}{3} = coder.typeof(ARGS_1_3);
ARGS{1}{4} = coder.typeof(ARGS_1_3);
ARGS{1}{5} = coder.typeof(0);
ARGS{1}{6} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg SimulateSchlogl2cell -args ARGS{1} -nargout 9

