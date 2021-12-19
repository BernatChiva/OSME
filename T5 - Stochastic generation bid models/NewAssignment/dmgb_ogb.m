function [] = dmgb_ogb(range_t)
% F.-Javier Heredia f.javier.heredia(at)upc.edu http://gnom.upc.edu/heredia
% Code under Creative Commons Attribution-NonCommercial-NoDerivs
% 3.0 Unported License http://creativecommons.org/licenses/by-nc-nd/3.0/
% 
% script to plot the optimal generation bid for model DMBG-R and DMBG-RB
%
% range_t: list of time periods to plot
%
lbgAcc = importdata('dmgb_lbg.a2m');
pbgAcc = importdata('dmgb_pbg.a2m');

pgb    = importdata('dmgb_pgb.a2m');

pbgAcc = pbgAcc';
lbgAcc = lbgAcc';

[nbg,tp] = size(lbgAcc);
[lbgAcc_s,IXg] = sort(lbgAcc);

for t = 1:tp
    pbgAcc_s(:,t) = pbgAcc(IXg(:,t),t); 
end
xg = zeros(nbg+2,tp); yg = zeros(nbg+2,tp);

for t = 1:tp
    for j = 2:nbg+1
        yg(j,t) = lbgAcc_s(j-1,t);
        if j==2
            xg(j,t) = 0;
        else
            xg(j,t) = xg(j-1,t) + pbgAcc_s(j-2,t);
        end
    end
    xg(nbg+2,t) = xg(nbg+1,t) + pbgAcc_s(nbg,t);
    yg(nbg+2,t) = yg(nbg+1,t);
end

tu = importdata('dmgb_tu.a2m');
cq=tu(1);
cl=tu(2);
pgmin = tu(3);
pgmax = tu(4);
laDmax = tu(5);
pgmin_mprice = 2*cq*pgmin + cl;
pgmax_mprice = 2*cq*pgmax + cl;

energy = max(xg(:));
price  = max(yg(:));
maxenergy = 1.05*pgmax;
maxprice = 1.05*price;

ntp = size(range_t,2);
nrows = 2;
ncol = ceil(ntp/nrows);
if ntp <= ncol nrows = 1; end
nsp = 0;
for t=range_t
    nsp = nsp+1;
%    subplot(ceil(tp/ncol),ncol,nsp)
    subplot(nrows,ncol,nsp)
    x = zeros(nbg,1);
    y = zeros(nbg,1);
    x(1) = pbgAcc_s(1,t);
    y    = lbgAcc_s(:,t);
    for i=2:nbg
        x(i) = x(i-1) + pbgAcc_s(i,t);
    end
    hold on;
    box on;
    xlabel({'{\itP} [MW]',['\fontname{arial}{\it\bf t = }',num2str(t)]},'Fontname','Cambria Math','Fontsize',20)
%    title({['\fontname{arial}{\it\bf t = }',num2str(t)]},'FontName','Cambria Math','Fontsize',20,'FontWeight','bold','HorizontalAlignment','left');
%    ax = gca;ax.TitleHorizontalAlignment = 'left';
    if ceil((t-1)/ncol) == (t-1)/ncol
        ylabel('{\it\lambda^D} [€/MWh]','Fontname','Cambria Math','Fontsize',20)
    end
    if t > (tp/ncol-1)*ncol
        %xlabel('{\itP} [MW]','Fontname','Cambria Math','Fontsize',20)
        xlabel({'{\itP} [MW]',['\fontname{arial}{\it\bf t = }',num2str(t)]},'Fontname','Cambria Math','Fontsize',20)
    end
    set(gcf,'Color','w'); set(gca,'Fontsize',20);
    if maxprice > 0
        axis ([0.01 maxenergy 0 maxprice]);
    else
        axis ([0.01 maxenergy 0 laDmax]);
    end
    energy = max(xg(1:nbg+2,t));
    price  = max(yg(1:nbg+2,t));
    gran = 100;
    d = 10^-4;
    if energy > 0
        x_pgmax_mprice = 0:maxenergy/gran:maxenergy;
        y_pgmax_mprice = pgmax_mprice:d/gran:pgmax_mprice+d;
        x_pgmax_pgmax  = pgmax:d/gran:pgmax+d;
        y_pgmax_pgmax  = 0:maxprice/gran:maxprice;

        x_pgmin_mprice = 0:maxenergy/gran:maxenergy;
        y_pgmin_mprice = pgmin_mprice:d/gran:pgmin_mprice+d;
        x_pgmin_pgmin  = pgmin:d/gran:pgmin+d;
        y_pgmin_pgmin  = 0:maxprice/gran:maxprice;

        x_laDmax_la = 0:maxenergy/gran:maxenergy;
        y_laDmax_la = laDmax:d/gran:laDmax+d;
        x_laDmax_pgmax = pgmax;
%        y_laDmax_pgmax = 0:laDmax/gran:laDmax;
        y_laDmax_pgmax = 0:maxprice/gran:maxprice;

        xmarginal = pgmin:(pgmax-pgmin)/gran:pgmax;
        ymarginal = 2*cq*xmarginal+cl;
        
        plot(xmarginal,ymarginal,'-k','Color',[.23,.54,.57],'LineWidth',1)
        plot(x_pgmax_mprice,y_pgmax_mprice,':k','LineWidth',1)
        plot(x_pgmin_mprice,y_pgmin_mprice,':k','LineWidth',1)
        plot(x_pgmin_pgmin,y_pgmin_pgmin,':k','LineWidth',1)
        plot(x_pgmax_pgmax,y_pgmax_pgmax,':k','LineWidth',1)
        % Pmin, Pmax
        text(pgmin, maxprice+7,'\itP^{min}_G','Horizontalalignment','center','Fontsize',16)
        text(pgmax, maxprice+7,'\itP^{max}_G','Horizontalalignment','center','Fontsize',16)
        if laDmax >0
            plot(x_laDmax_pgmax,y_laDmax_pgmax,':k','LineWidth',6)
            plot(x_laDmax_la,y_laDmax_la,':k','LineWidth',1)
            if ceil((nsp)/ncol) == (nsp)/ncol || nsp == ntp
                text(maxenergy+5,laDmax,'\it\lambda^{max}','Horizontalalignment','left','Fontsize',16)
            end
        else
            plot(x_pgmax_pgmax,y_pgmax_pgmax,':k','LineWidth',6)
        end
        if ceil((nsp)/ncol) == (nsp)/ncol || nsp == ntp
            text(maxenergy+5,pgmax_mprice,'\it\lambda^m_G(P_G^{max})','Horizontalalignment','left','Fontsize',16)
            text(maxenergy+5,pgmin_mprice,'\it\lambda^m_G(P_G^{min})','Horizontalalignment','left','Fontsize',16)
        end
    end
    
    % Box of the OGB
    if max(pgb > 0)
        gran1 = 100;
        d=10^-4;
        xl = pgb(t):d/gran1:pgb(t)+d; yl = maxprice*(xl-pgb(t))/d;%0:maxprice/gran1:maxprice;
        xr = (pgmax-d):d/gran1:pgmax; yr = maxprice*(xr-(pgmax-d))/d;%0:maxprice/gran1:maxprice;
        xu = pgb(t):pgmax/gran1:pgmax; yu = -d*(xu-pgb(t))+maxprice;
        xb = xu; yb = d*xb;
        plot(xl,yl,'--','Color',[.17,.17,0.53],'LineWidth',2);
        plot(xr,yr,'--','Color',[.17,.17,0.53],'LineWidth',2);
        if energy >0
            plot(xu,yu,'--','Color',[.17,.17,0.53],'LineWidth',2);
            plot(xb,yb,'--','Color',[.17,.17,0.53],'LineWidth',2);
        end
    end
    % OGB
    if ntp <= 4
        fsogb = 6;
    else
        fsogb = 4;
    end
    xg(:,t) = xg(:,t) + pgb(t);
    stairs(xg(:,t),yg(:,t),'-','Color',[.17,.17,0.53],'LineWidth',fsogb)

    % Scenarios
    x(:) = x(:) + pgb(t);
    plot(x(2:nbg-1),y(2:nbg-1),'mo','Color',[.8,.0,0.8],'MarkerSize',fsogb-1,'LineWidth',fsogb-1);
    % PGB
    xx=0:1:pgb(t); yy=xx*10^-10;
    plot(xx,yy,'LineWidth',fsogb, 'Color',[0,0.75,0]);
    plot(pgb(t),0,'mo','MarkerSize',fsogb,'LineWidth',fsogb, 'Color',[0,0.75,0]);
    if pgb(t)>0
        text(pgb(t)/2, 5,'\it P^{B*}_G','Horizontalalignment','center','Fontsize',16, 'Color',[0,0.75,0])
    end
    plot(pgb(t),0,'mo','MarkerSize',fsogb,'LineWidth',fsogb, 'Color',[0,0.75,0]);
    %
    hold off;
end
%subplot (ceil(tp/ncol),ncol,tp-1);
% text(maxenergy/2, -8, ['\bf{\it SW ^*} = ',num2str(SWT)],'Fontname','Cambria Math','Fontsize',20);

end

