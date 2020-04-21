function ScalingTwocellSchloglCommandline_diff_params(savedir, theta_x, theta_y, g_x, g_y, h_x, h_y, log10nc_x, log10nc_y)
% Simulates linear birth-death model to test simulation

rng('shuffle'); % Essential !!! Shuffle random number generator

getTranjectories = false; % True just get trajectories, false gets distributions
% Thing is, for large nc, the joint distribution is too large for mex !


nc_x = round(10^log10nc_x);
nc_y = round(10^log10nc_y);

mkdir(savedir);


% Choose Ising variables,
disp(['Theta_x ' num2str(theta_x) ' ; Theta_y ' num2str(theta_y) ...
      ' ; h_x ' num2str(h_x) ' ; h_y ' num2str(h_y) ...
      ' ; nc ' num2str(nc_x)]);
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
tspan = [1E5, 2E6];  % Actual time
n0 = [Ising_x.nc, Ising_y.nc];    % Initial copy number

if(~getTranjectories)
    savefile = [savedir filesep '/out__nc_' num2str(nc_x) '__thetax_' num2str(Ising_x.theta) '__thetay_' num2str(Ising_y.theta) '__g_' num2str(Ising_x.g) ...
        '__hx_' num2str(h_x) '__hy_' num2str(h_y) '.mat'];
else
  %    savefile = [savedir filesep '/Trajectories_out_nc_' num2str(nc_x) '_theta_' num2str(Ising_x.theta) '_g_' num2str(Ising_x.g) '.mat'];
end

if isfile(savefile)
    disp(['File already exists. Abortin. ' savefile]);
    return 
end

Schlogl_x = SchloglFromIsing(Ising_x);
Schlogl_y = SchloglFromIsing(Ising_y);

GillespieOut = struct;
GillespieOut.getTranjectories = getTranjectories;
GillespieOut.Ising_x = Ising_x;
GillespieOut.Ising_y = Ising_y;
GillespieOut.Schlogl_x = Schlogl_x;
GillespieOut.Schlogl_y = Schlogl_y;


if(~getTranjectories)
    if(Ising_x.h + Ising_y.h < -0.6)
        disp(['Ignoring h_x+h_y=' num2str(Ising_x.h + Ising_y.h) ' < -0.6'])
        return
    end

    [nSteps,Pn,Pm,Pnm,Pz,Pw,tau_n, tau_m, batchMeans] = ...
        SimulateSchlogl2cell_mex(tspan, n0, Schlogl_x, Schlogl_y, 5E8, 1E7);
    GillespieOut.Pn = Pn;
    GillespieOut.Pm = Pm;
    GillespieOut.Pnm = sparse(Pnm);
    GillespieOut.Pz = Pz;
    GillespieOut.Pw = Pw;
else % Trajectories
    [nSteps,T,X,tau_n, tau_m, batchMeans] = ...
        SimulateSchlogl2cellTrajectories_mex(tspan, n0, Schlogl_x, Schlogl_y, 1E8, 1E6);
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
