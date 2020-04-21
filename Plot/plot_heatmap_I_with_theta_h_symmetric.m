% plot I(\theta,h) for symmetric cells
% k_{1}^{1} = 1 for both cells
% g = 1 for both cells
% theta & h are varied for both cells, but are always equal,
% i.e. theta1 = theta2, h1 = h2 for every data point
% The color is the Shannon mutual information, in units of nats
% The bottom left corner is cutoff because of the restriction
% h + theta > 1/3

fulltab = readtable('../Data/symmetric_h_theta_nc3000-collected.csv','delimiter',',');
fulltab = fulltab(fulltab.nc_x==3000, :);

symmtab = fulltab(fulltab.h_x == fulltab.h_y,:);
symmtab = symmtab(symmtab.theta_x == symmtab.theta_y, :);

h_symm = unique(symmtab.h_x);
theta_symm = unique(symmtab.theta_x);

info_symm = zeros(length(theta_symm), length(h_symm));
for tt=1:length(theta_symm)
    for hh=1:length(h_symm)
        g = symmtab(symmtab.h_x == h_symm(hh) & symmtab.theta_x == theta_symm(tt), :);
%         if height(g) ~= 1
%             disp(['Warning! g has ' num2str(height(g)) ' rows !']);
%             continue
%         end
        if theta_symm(tt)<-0.15
            info_symm(tt, hh) =0;
        else
            info_symm(tt, hh) = g.I;
        end
    end
end

newfigure(4.5,3.75);
% surfc(h_symm,theta_symm,info_symm);
% imagesc(h_symm,theta_symm,log10(info_symm));
imagesc(h_symm,theta_symm,info_symm);
colormap(hot(256));
colorbar;
set(gca,'YDir','Normal')
xlabel('$h$','Interpreter','latex');
ylabel('$\theta$','Interpreter','latex');
set(gca,'FontSize',18);
xlim([-0.1, 0.1]);
ylim([-0.1, 0.1]);
t = text(0.127, 0.107, '$I$', 'Interpreter', 'Latex', 'FontSize', 18);
print(gcf,'-dpng','ShannonI_heatmap_theta_h.png','-r600');


