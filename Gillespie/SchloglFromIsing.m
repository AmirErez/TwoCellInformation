function Schlogl = SchloglFromIsing(Ising)
   Schlogl = struct;
   Schlogl.s = 3*(Ising.nc-1);
   Schlogl.K2 = 3*Ising.theta*Ising.nc^2 + 3*(Ising.nc-1)*(Ising.nc-2)+1;
   Schlogl.K = sqrt(Schlogl.K2);
   Schlogl.a = ((3*Ising.theta+3*Ising.h+1)*Ising.nc^3 - 6*(Ising.nc-1))/(Schlogl.K^2);
   Schlogl.N = Ising.nc*10;
   if(isfield(Ising, 'g'))
       Schlogl.g = Ising.g;
   end
