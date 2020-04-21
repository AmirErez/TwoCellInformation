% Used to create Fig.2 in the manuscript

figname = 'Fig3ab';
newfigure(3.375/2,3.375/1.8);
set(gca,'FontSize',9);
% nmap = colormap(parula(length(collected.ncs)+5));
% [ha, pos] = tight_subplot(2,2,[.05 .08],[.08 .01],[.08 .01]);

% manythetasfile = '../Data/ProductionManyThetas.csv';
many_thetas_file = '../Data/ProductionManyThetasParticularNcs.csv';
theta_zero_file = '../Data/ProductionTheta0Scaling.csv';

tabManyThetas = readtable(many_thetas_file);
tabManyThetas = sortrows(tabManyThetas,'nc');
tabThetaZero = readtable(theta_zero_file);
tabThetaZero = sortrows(tabThetaZero,'nc');
tabThetaZero.g = round(log10(tabThetaZero.g)*10)/10;

% (a) Shannon entropy unscaled
subplot(2,1,1);
set(gca,'FontSize',9);
hold on
subtab=tabManyThetas(tabManyThetas.g==1,:);
ncs = unique(tabManyThetas.nc);
nmap = colormap(parula(length(ncs)+3));

text(-0.18,1.75,'(a)');

for nn=1:length(ncs)
   tabnc = subtab(subtab.nc==ncs(nn),:); 
   tabnc = sortrows(tabnc,'theta');
   f = find(tabnc.IShannon);   
   plot(tabnc.theta(f),tabnc.IShannon(f),'-','LineWidth', 1,'Color',nmap(nn,:),'DisplayName',['n_c = ' num2str(ncs(nn))]);   
end
xlabel('$\theta$','Interpreter','latex');
ylabel('$I$','Interpreter','latex');

plot([0,0],[0,2],'--k');
% text(0.55,-0.05,'$\theta$','Interpreter','latex');
% text(-0.18,2,'$I_2$','Interpreter','latex');
c = colorbar();
c.TickLabels={'$10^{2.5}$', '$n_c$', '$10^{3.7}$'};
c.TickLabelInterpreter = 'latex';
mxy = max(subtab.IShannon);
ylim([0,mxy]);
xlim([-0.2,0.4]);
set(gca,'XTick',[-0.2,0,0.2,0.4]);
% set(gca,'XTickLabel',[0,0.4]);
set(gca,'YTick',[0:1:2]);
set(gca,'YTickLabel',[0:1:2]);

% axes('Position',[0.55 0.65 0.35 0.32])
% box on
% set(gca,'FontSize',9);
% hold on
% for nn=1:length(ncs)
%    tabnc = subtab(subtab.nc==ncs(nn),:); 
%    tabnc = sortrows(tabnc,'theta');
%    f = find(tabnc.IShannon);   
%    plot(tabnc.theta(f),tabnc.IShannon(f),'.-','Color',nmap(nn,:),'DisplayName',['n_c = ' num2str(ncs(nn))]);   
% end
% plot([0,0],[0.45,mxy+0.2],'--k','LineWidth',1);
% xlabel('$\theta$','Interpreter','latex');
% ylabel('$I$','Interpreter','latex');
% xlim([-0.05, 0.05]);
% ylim([0.45, mxy+0.2]);
% set(gca,'YTickLabel',[]);
% set(gca,'XTick',[-0.05, 0.05]);

%% Fig 3b: Show of I(theta=0, g) with nc

subplot(2,1,2);
set(gca,'FontSize',9);
hold on

% gs = unique(tabThetaZero.g);
gs = [-1:0.5:0.5];
ncs = unique(tabManyThetas.nc);
gmap = colormap(summer(length(gs)+1));

for gg=length(gs):-1:1
    if(gs(gg)==-inf), continue; end
    subtab = tabThetaZero(tabThetaZero.g==gs(gg),:);
    subtab = subtab(subtab.theta==0,:);
    subtab = sortrows(subtab,'nc');
    plot(subtab.nc,subtab.IShannon,'.-','Color',gmap(gg,:),'DisplayName',[num2str(gs(gg))]);
end
[lb,licons,lplots,ltxt]= legend('show');
lb.Box = 'Off';
lb.Position = [0.72        0.33       0.1      0.10082];
% licons(1).Position(1) = 0.6;
% licons(2).Position(1) = 0.5

semilogx(10.^[2,4], -0.25+0.25*log(10.^[2,4]),'--k');

% set(gca,'YScale','log');
xticks(10.^[2,3,4]);
set(gca,'XScale','log');
xlabel('$n_c$','Interpreter','Latex');
ylabel('$I(n_c, g ; \theta=0)$','Interpreter','Latex');
xlim(10.^[2,5]);
ylim([0,2]);
text(10^2.1,2.2,'(b)');
print(gcf,'-dpng','Fig3ab.png','-r600');


