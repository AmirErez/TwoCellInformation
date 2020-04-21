function Hill = HillFromIsing(Ising, H)
   Hill = struct;
   Hill.a = Ising.nc*(H-1)*((H+1)*(H+1)*(Ising.theta+Ising.h)+4)/(H+1)/((H*H-1)*Ising.theta +4);
   Hill.s = Ising.nc*16*H/(H*H-1)/((H*H-1)*Ising.theta+4);
   Hill.K = Ising.nc*((H+1)/(H-1))^(1/H);
   Hill.H = H
   Hill.N = Ising.nc*10;
   if(isfield(Ising, 'g'))
       Hill.gamma = Ising.g*(H*H-1)/(Ising.theta*(H*H-1)+4);
   end
