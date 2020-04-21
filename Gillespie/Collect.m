% Collects Gillespie output in folder and makes table with params and corr time
% all_collect = {'symmetric_h_theta_nc3000'};
% all_collect = {'scan_thetax_thetay_h0_nc3000-hill'};
% all_collect = {'scan_thetax_thetay_h0_nc3000-schlogl'};
% all_collect = {'scan_hx_hy_theta0_nc3000-hill'};
% all_collect = {'scan_hx_hy_theta0_nc3000'};
all_collect = {'randparams0_1-nc_3000-schlogl'};

for aa=1:length(all_collect)
    disp(['Doing ' all_collect{aa} ]);

    collectdir = all_collect{aa};

    d = dir([collectdir filesep '*.mat']);
    
    tabInfo = table;
    
    for ff=1:length(d)
        disp(d(ff).name);
        fullname = [collectdir filesep d(ff).name];
        try
           load(fullname);
           loaded = GillespieOut; 
           tempstruct = struct;
           tempstruct.theta_x = loaded.Ising_x.theta;
           tempstruct.theta_y = loaded.Ising_y.theta;
           tempstruct.h_x = loaded.Ising_x.h;
           tempstruct.h_y = loaded.Ising_y.h;
           tempstruct.nc_x = loaded.Ising_x.nc;
           tempstruct.nc_y = loaded.Ising_y.nc;
           tempstruct.n_steps = loaded.nSteps;
           tempstruct.tau_x = loaded.tau_n;
           tempstruct.tau_y = loaded.tau_m;
%           tempstruct.filename = fullname;
        catch
           warning(['Failed loading ' fullname ]);
           continue
        end
        Pn = loaded.Pn/sum(loaded.Pn);
        tempstruct.mean_n = sum((0:length(Pn)-1)'.*Pn);
        Pm = loaded.Pm/sum(loaded.Pm);
        tempstruct.mean_m = sum((0:length(Pn)-1)'.*Pn);
%        Pxy = loaded.Pnm/sum(sum(loaded.Pnm));
        v = nonzeros(loaded.Pn);
        v = v / sum(v);
        tempstruct.S_x = -sum(v.*log(v));
        v = nonzeros(loaded.Pm);
        v = v / sum(v);
        tempstruct.S_y = -sum(v.*log(v));
        v = nonzeros(loaded.Pnm);
        v = v / sum(v);
        tempstruct.S_xy = -sum(v.*log(v));

        margx = nonzeros(sum(loaded.Pnm, 1));
        margx = margx / sum(margx);
        margy = nonzeros(sum(loaded.Pnm, 2));
        margy = margy / sum(margy);

        tempstruct.I_try = tempstruct.S_x + tempstruct.S_y - tempstruct.S_xy;
        tempstruct.I = -sum(margx.*log(margx))-sum(margy.*log(margy))-tempstruct.S_xy;
        if(abs(imag(tempstruct.I))>1E-6)
           tempstruct.I = 0;
        else
           tempstruct.I = real(tempstruct.I);
        end
        tabInfo = [tabInfo; struct2table(tempstruct)];
    end
    disp('Finished');
    writetable(tabInfo,[collectdir filesep 'collected.csv']);
end
