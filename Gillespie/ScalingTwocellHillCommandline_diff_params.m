function ScalingTwocellHillCommandline_diff_params(savedir, theta_x, theta_y, g_x, g_y, h_x, h_y, log10nc_x, log10nc_y, H)
% Simulates linear birth-death model to test simulation

rng('shuffle'); % Essential !!! Shuffle random number generator

getTranjectories = false; % True just get trajectories, false gets distributions
% Thing is, for large nc, the joint distribution is too large for mex !


nc_x = round(10^log10nc_x);
nc_y = round(10^log10nc_y);

mkdir(savedir);


% Choose Ising variables,
disp(['Theta ' num2str(theta_x) ' ; nc ' num2str(nc_x)]);
Ising_x.nc = nc_x;
Ising_x.theta = theta_x;
Ising_x.h = h_x;
Ising_x.g = g_x;
% Ising_x.tspan = [0, 2E5];  % Actual time
Ising_x.tspan = [0, 0];  % IGNORE Actual time

Ising_y.nc = nc_y;
Ising_y.theta = theta_y;
Ising_y.h = h_y;
Ising_y.g = g_y;

% tspan = [0, 2E5];  % Actual time
% tspan = [1E5, 1E6];  % Actual time
tspan = [1E5, 5E5];  % Actual time
n0 = [Ising_x.nc, Ising_y.nc];    % Initial copy number

if(~getTranjectories)
    savefile = [savedir filesep '/out__nc_' num2str(nc_x) '__thetax_' num2str(Ising_x.theta) '__thetay_' num2str(Ising_y.theta) '__g_' num2str(Ising_x.g) ...
        '__hx_' num2str(h_x) '__hy_' num2str(h_y) '.mat'];
else
    savefile = [savedir filesep '/Trajectories_out_nc_' num2str(nc_x) '_theta_' num2str(Ising_x.theta) '_g_' num2str(Ising_x.g) '.mat'];
end

if isfile(savefile)
    disp(['File already exists. Abortin. ' savefile]);
    return 
end

Hill_x = HillFromIsing(Ising_x, H);
Hill_y = HillFromIsing(Ising_y, H);

GillespieOut = struct;
GillespieOut.getTranjectories = getTranjectories;
GillespieOut.Ising_x = Ising_x;
GillespieOut.Ising_y = Ising_y;
GillespieOut.Hill_x = Hill_x;
GillespieOut.Hill_y = Hill_y;


if(~getTranjectories)
    [nSteps,Pn,Pm,Pnm,~,~,tau_n, tau_m, batchMeans] = ...
        SimulateHill2cell_mex(tspan, n0, Hill_x, Hill_y, 5E8, 1E7);
    GillespieOut.Pn = Pn;
    GillespieOut.Pm = Pm;
    GillespieOut.Pnm = sparse(Pnm);
    GillespieOut.Pz = [];
    GillespieOut.Pw = [];
else % Trajectories
    [nSteps,T,X,tau_n, tau_m, batchMeans] = ...
        SimulateSchlogl2cellTrajectories_mex(tspan, n0, Hill_x, Hill_y, 1E8, 1E6);
    GillespieOut.T = T;
    GillespieOut.X = int16(X);
end


GillespieOut.nSteps = nSteps;
GillespieOut.tau_n = tau_n;
GillespieOut.tau_m = tau_m;
GillespieOut.batchMeans = batchMeans;

save(savefile,'GillespieOut');
clear('GillespieOut');

disp('Finished successfully');
