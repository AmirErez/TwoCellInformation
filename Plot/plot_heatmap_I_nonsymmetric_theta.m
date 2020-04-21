% plot I(\theta,h) for symmetric cells
% k_{1}^{1} = 1 for both cells
% g = 1 for both cells
% theta & h are varied for both cells, but are always equal,
% i.e. theta1 = theta2, h1 = h2 for every data point
% The color is the Shannon mutual information, in units of nats
% The bottom left corner is cutoff because of the restriction
% h + theta > 1/3


fulltab = readtable('../Data/scan_thetax_thetay_h0_nc3000-schlogl-collected.csv','delimiter',',');
fulltab = fulltab(fulltab.nc_x==3000, :);
hzerotab = fulltab(fulltab.h_x == 0 & fulltab.h_y == 0,:);

theta1_varytheta = unique(hzerotab.theta_x);
theta2_varytheta = unique(hzerotab.theta_y);

info_varytheta = zeros(length(theta2_varytheta), length(theta1_varytheta));
for t2=1:length(theta2_varytheta)
    for t1=1:length(theta1_varytheta)
        g = hzerotab(hzerotab.theta_x == theta1_varytheta(t1) & hzerotab.theta_y == theta2_varytheta(t2), :);
        if height(g) == 0
            continue
        end
        if height(g) ~= 1
            disp(['Warning! g has ' num2str(height(g)) ' rows !']);
            continue
        end
        info_varytheta(t2, t1) = g.I;
    end
end

newfigure(4.5,3.75);
% surfc(theta2_varytheta,theta1_varytheta,info_varytheta);
% imagesc(theta2_varytheta,theta1_varytheta,log10(info_varytheta));
imagesc(theta1_varytheta,theta2_varytheta,info_varytheta);
colormap(hot(256));
colorbar;
% view(2);
set(gca,'YDir','Normal');
xlabel('$\theta_X$','Interpreter','Latex');
ylabel('$\theta_Y$','Interpreter','Latex');
set(gca,'FontSize',18);

hold on
warning('Setting g=1');
g=1;
plot(theta1_varytheta, -g*theta1_varytheta./(theta1_varytheta+g), '--b', 'LineWidth', 2);
xticks([-0.3, 0, 0.3]);
yticks([-0.3, 0, 0.3]);
t = text(0.38, 0.33, '$I$', 'Interpreter', 'Latex', 'FontSize', 18);
t1 = text(-0.52, 0.32, '(a)', 'Interpreter', 'Latex', 'FontSize', 22);
print(gcf,'-dpng','ShannonI_nonsymmetric_theta.png','-r600');

%% Hill dynamics with constant gamma not g

fulltab = readtable('../Data/scan_thetax_thetay_h0_nc3000-hill-collected.csv','delimiter',',');
fulltab = fulltab(fulltab.nc_x==3000, :);
hzerotab = fulltab(fulltab.h_x == 0 & fulltab.h_y == 0,:);

theta1_varytheta = unique(hzerotab.theta_x);
theta2_varytheta = unique(hzerotab.theta_y);

info_varytheta = zeros(length(theta2_varytheta), length(theta1_varytheta));
for t2=1:length(theta2_varytheta)
    for t1=1:length(theta1_varytheta)
        g = hzerotab(hzerotab.theta_x == theta1_varytheta(t1) & hzerotab.theta_y == theta2_varytheta(t2), :);
        if height(g) ~= 1
            disp(['Warning! g has ' num2str(height(g)) ' rows !']);
            continue
        end
        info_varytheta(t2, t1) = g.I;
        if (theta2_varytheta(t2)<=-0.25 && theta1_varytheta(t1)<=-0.25)
            info_varytheta(t2, t1) = 0;
        end
    end
end

newfigure(4.5,3.75);
% surfc(theta2_varytheta,theta1_varytheta,info_varytheta);
% imagesc(theta2_varytheta,theta1_varytheta,log10(info_varytheta));
imagesc(theta1_varytheta,theta2_varytheta,info_varytheta);
colormap(hot(256));
colorbar;
% view(2);
set(gca,'YDir','Normal');
xlabel('$\theta_X$','Interpreter','Latex');
ylabel('$\theta_Y$','Interpreter','Latex');
set(gca,'FontSize',18);

hold on
warning('Setting g=1');
g=1;
plot(theta1_varytheta, -g*theta1_varytheta./(theta1_varytheta+g), '--b', 'LineWidth', 2);
xlim([-0.3, 0.3]);
ylim([-0.3, 0.3]);
xticks([-0.3, 0, 0.3]);
yticks([-0.3, 0, 0.3]);
t = text(0.38, 0.33, '$I$', 'Interpreter', 'Latex', 'FontSize', 18);
t1 = text(-0.52, 0.32, '(a)', 'Interpreter', 'Latex', 'FontSize', 22);

print(gcf,'-dpng','ShannonI_nonsymmetric_theta_Hill.png','-r600');