% SIMULATEHILL2CELL_CODEGENSCRIPT   Generate MEX-function SimulateHill2cell_mex
%  from SimulateHill2cell.
% 
% Script generated from project 'SimulateHill2cell.prj' on 03-Jan-2020.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

%% Define argument types for entry-point 'SimulateHill2cell'.
ARGS = cell(1,1);
ARGS{1} = cell(6,1);
ARGS{1}{1} = coder.typeof(0,[1 2]);
ARGS{1}{2} = coder.typeof(0,[1 2]);
ARGS_1_3 = struct;
ARGS_1_3.a = coder.typeof(0);
ARGS_1_3.s = coder.typeof(0);
ARGS_1_3.K = coder.typeof(0);
ARGS_1_3.H = coder.typeof(0);
ARGS_1_3.N = coder.typeof(0);
ARGS_1_3.gamma = coder.typeof(0);
ARGS{1}{3} = coder.typeof(ARGS_1_3);
ARGS_1_4 = struct;
ARGS_1_4.a = coder.typeof(0);
ARGS_1_4.s = coder.typeof(0);
ARGS_1_4.K = coder.typeof(0);
ARGS_1_4.H = coder.typeof(0);
ARGS_1_4.N = coder.typeof(0);
ARGS_1_4.gamma = coder.typeof(0);
ARGS{1}{4} = coder.typeof(ARGS_1_4);
ARGS{1}{5} = coder.typeof(0);
ARGS{1}{6} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg SimulateHill2cell -args ARGS{1} -nargout 9

