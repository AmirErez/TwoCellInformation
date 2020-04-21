function [nSteps,Pn,Pm,Pnm,Pz,Pw,tau_n, tau_m, batchMeans] = SimulateSchlogl2cell(tspan, x0,...
                                  Schlogl_n, Schlogl_m, MAX_OUTPUT_LENGTH,DISPLAY_EVERY)
%#codegen
%DIRECTMETHOD Implementation of the Direct Method variant of the Gillespie algorithm
%   Returns:
%       t:              time vector          (Nreaction_events x 1)
%       x:              species amounts      (Nreaction_events x Nspecies)    
%
%   Required:
%       tspan:          Initial and final times, [t_init, t_final].
%
%       x0:             Initial species amounts, [S1_0, S2_0, ... ].
%
%       stoich_matrix:  Matrix of stoichiometries (Nreactions x Nspecies).
%                       Each row gives the stoichiometry of a reaction.
%
%
%   Reference: 
%       Gillespie, D.T. (1977) Exact Stochastic Simulation of Coupled
%       Chemical Reactions. J Phys Chem, 81:25, 2340-2361.
%
%   Adapted by Amir Erez 2017-2020 (amir.b.erez@gmail.com)

%% Initialize
%num_rxns = size(stoich_matrix, 1);
isFullTimeseries = false;
% if(isFullTimeseries)
%     T = zeros(MAX_OUTPUT_LENGTH, 1);
%     X = zeros(MAX_OUTPUT_LENGTH, 2);
%     T(1)     = tspan(1);
%     X(1,:)   = x0;
% else
%    t = [];
%    x = [];
% end

prevX = x0;
newX = NaN*prevX;
curT = 0;
rxn_count = 1;

% For correlation time using method of batch means
%nBatches = floor(tspan(2)^(1/3));
%batchTime = tspan(2)/nBatches;
batchTime = 1000;
nBatches = floor(tspan(2)/batchTime);

batchMeans = zeros(nBatches,2);
totalVariance = zeros(1,2);
running_batch_sum = zeros(1,2);
running_batch_sq_sum = zeros(1,2);
% batchVals = NaN(5E8,2);

% DISPLAY_EVERY = 100000;

fprintf('Running with params\n');
fprintf('n: = a=%.2f, g=%.2f, K=%.2f, K2=%.2f, N=%f, s=%.2f\n',Schlogl_n.a, Schlogl_n.g, Schlogl_n.K, Schlogl_n.K2, Schlogl_n.N, Schlogl_n.s);
fprintf('m: = a=%.2f, g=%.2f, K=%.2f, K2=%.2f, N=%f, s=%.2f\n',Schlogl_m.a, Schlogl_m.g, Schlogl_m.K, Schlogl_m.K2, Schlogl_m.N, Schlogl_m.s);
fprintf('-------------------------------------------------------------------\n');

nc_n = Schlogl_n.s/3+1
nc_m = Schlogl_m.s/3+1

k_n1minus = 1;
k_m1minus = k_n1minus/(Schlogl_n.K^2)*(Schlogl_m.K^2);

k_n1plus = Schlogl_n.a*k_n1minus;
k_m1plus = Schlogl_m.a*k_m1minus;
k_n2minus = k_n1minus/(Schlogl_n.K^2);
k_m2minus = k_m1minus/(Schlogl_m.K^2);
k_n2plus = k_n2minus*Schlogl_n.s;
k_m2plus = k_m2minus*Schlogl_m.s;
gamma_mn = Schlogl_m.g*k_m1minus/3*((Schlogl_m.s+3)/Schlogl_m.K)^2;
gamma_nm = Schlogl_n.g*k_n1minus/3*((Schlogl_n.s+3)/Schlogl_n.K)^2;


Pn = zeros(x0(1)*10,1);
Pm = zeros(x0(1)*10,1);
Pnm = zeros(x0(1)*6,x0(1)*6);
Pz = zeros(x0(1)*20,1);
Pw = zeros(x0(1)*10,1);

stoich_matrix = [
1 0
-1 0
-1 1
1 -1
0 1
0 -1
];

%% MAIN LOOP
batchCounter = 0;
counterInBatch = 0;
is_burnin = true;
while curT < tspan(2)        

    if(curT > batchCounter*batchTime)
        if(batchCounter>0 & counterInBatch>0) % Statistics for previous batch
%            batchMeans(batchCounter,:) = nanmean(batchVals(1:counterInBatch,1));
%            totalVariance = (totalVariance*batchCounter + nanmean(batchVals(1:counterInBatch,1).^2) - batchMeans(batchCounter,:).^2)/(batchCounter+1);
             batchMeans(batchCounter,:) = running_batch_sum/counterInBatch;
             batch_variance = running_batch_sq_sum/counterInBatch - batchMeans(batchCounter,:).^2 ;
             totalVariance = (totalVariance*batchCounter + batch_variance)/(batchCounter+1);
             
        end        
        batchCounter = batchCounter+1;
        counterInBatch = 0;
        running_batch_sum = zeros(1,2);
        running_batch_sq_sum = zeros(1,2);
    end
    

    if(is_burnin==true && curT>=tspan(1))
          is_burnin = false;
          fprintf('Finished burnin\n');
          curT = 0;
          rxn_count = 1;
    end    
    if(mod(rxn_count,DISPLAY_EVERY)==0)
         if(is_burnin==true)
            fprintf('Burnin t = %.2f\n',curT);
         else 
            fprintf('t = %.2f\n',curT);
         end
    end

    % Calculate reaction propensities   
    a = [
        k_n1plus + k_n2plus*prevX(1)*(prevX(1)-1)
        k_n1minus*prevX(1) + k_n2minus*prevX(1)*(prevX(1)-1)*(prevX(1)-2)
        gamma_nm*prevX(1)
        gamma_mn*prevX(2)
        k_m1plus + k_m2plus*prevX(2)*(prevX(2)-1)
        k_m1minus*prevX(2) + k_m2minus*prevX(2)*(prevX(2)-1)*(prevX(2)-2)
        ];
    
    % Sample earliest time-to-fire (tau)
    a0 = sum(a);
    r = rand(1,2);
    tau = -log(r(1))/a0; %(1/a0)*log(1/r(1));
    
    % Sample identity of earliest reaction channel to fire (mu)
%     [~, mu] = histc(r(2)*a0, [0;cumsum(a(:))]); 
    
    % ...alternatively...
    mu = find((cumsum(a) >= r(2)*a0), 1,'first');
    
    % ...or...
    %mu=1; s=a(1); r0=r(2)*a0;
    %while s < r0
    %   mu = mu + 1;
    %   s = s + a(mu);
    %end
    
    if (isFullTimeseries && (rxn_count + 1 > MAX_OUTPUT_LENGTH))
        if(isFullTimeseries)
            t = T(1:rxn_count);
            x = X(1:rxn_count,:);
        end
        disp('Number of reaction events exceeded the number pre-allocated. Simulation terminated prematurely.');
        break;
    end
    
    % Update time and carry out reaction mu
    newX = prevX + stoich_matrix(mu,:); 
    curT = curT + tau;
    if(isFullTimeseries)
        T(rxn_count+1)   = curT;
        X(rxn_count+1,:) = newX;
    end

%    batchVals(counterInBatch,:) = newX;
    running_batch_sum = running_batch_sum + newX;
    running_batch_sq_sum = running_batch_sq_sum + newX.^2;
    counterInBatch = counterInBatch+1;
    
    Pn(prevX(1)+1) = Pn(prevX(1)+1)+tau;
    Pm(prevX(2)+1) = Pn(prevX(2)+1)+tau;
    Pnm(prevX(1)+1,prevX(2)+1) = Pnm(prevX(1)+1,prevX(2)+1)+tau;
    
    Pz(sum(prevX)+1) = Pz(sum(prevX)+1) + tau;
    w = abs(diff(prevX))+1;
    Pw(w) = Pw(w) + tau;
    
    rxn_count = rxn_count + 1;
    prevX = newX;
    
end  

Pn = Pn/sum(Pn);
Pm = Pm/sum(Pm);
Pnm = Pnm/sum(Pnm(:));
Pz = Pz / sum(Pz);
Pw = Pw / sum(Pw);
nSteps = rxn_count;

% Batch means:
% batchMeans(batchCounter,:) = nanmean(batchVals,1);
% totalVariance = (totalVariance*batchCounter + nanmean(batchVals.^2,1) - batchMeans(batchCounter,:).^2)/(batchCounter+1);
batchMeans = batchMeans(1:(batchCounter-1),:)
varBatches = mean(batchMeans.^2,1) - mean(batchMeans,1).^2;
tau_n = batchTime * varBatches(1)/2/totalVariance(1);
tau_m = batchTime * varBatches(2)/2/totalVariance(2);

% Return simulation time course
% t = T(1:rxn_count);
% x = X(1:rxn_count,:);
% if t(end) > tspan(2)
%     t(end) = tspan(2);
%     x(end,:) = X(rxn_count-1,:);
% end    

end

