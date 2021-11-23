function tcmpa_abf = tcmpa_abf(kS,tpS)
%
% F.-Javier Heredia f.javier.heredia(at)upc.edu http://gnom.upc.edu/heredia
% Code under Creative Commons Attribution-NonCommercial-NoDerivs
% 3.0 Unported License http://creativecommons.org/licenses/by-nc-nd/3.0/
% 
% Function to represent the NODAL CLEARING for model TCMPA for busses 
% in the list kS and time periods in the list tpS.
% The input data files are created by the AMPL script tcma2m.run
%
%
delimiterIn = ' ';
busL   = importdata("tcmpa_busL.a2m");
nB = size(busL,1);
nk  = size(kS,2);
kS  = sort(kS);
if nk == 0
    kS=[1:nB];
    nk = nB;
end
ntp = size(tpS,2);
tpS = sort(tpS);
tiledlayout(nk,ntp);
for k = kS
    lbgAcc = importdata(strcat(busL{k},'\tcmpa_lbG.a2m'));
    pbgAcc = importdata(strcat(busL{k},'\tcmpa_pbG.a2m'));
    PGAcc  = importdata(strcat(busL{k},'\tcmpa_PG.a2m'));
    bus    = importdata(strcat(busL{k},'\tcmpa_bus.a2m'));

    pbgAcc = pbgAcc';
    lbgAcc = lbgAcc';
    PGAcc = PGAcc';

    [nbg,tp] = size(lbgAcc);
    if ntp == 0
        tpS = [1:tp];
        ntp  = tp;
    end
    [lbgAcc_s,IXg] = sort(lbgAcc,1);

    clearvars pbgAcc_s PGAcc_s pbdAcc_s PDAcc_s;
    for t = 1:tp
        pbgAcc_s(:,t) = pbgAcc(IXg(:,t),t);
        PGAcc_s(:,t) = PGAcc(IXg(:,t),t);
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

    xmg = zeros(nbg+2,tp); ymg = zeros(nbg+2,tp);
    for t = 1:tp
        pmgAcc_s = PGAcc_s(PGAcc_s(:,t)>0,t);
        lmgAcc_s = lbgAcc_s(PGAcc_s(:,t)>0,t);
        [nmg,tpm] = size(lmgAcc_s);
        for j = 2:nmg+1
            ymg(j,t) = lmgAcc_s(j-1);
            if j==2
                xmg(j,t) = 0;
            else
                xmg(j,t) = xmg(j-1,t) + pmgAcc_s(j-2);
            end
        end
        if nmg > 0
            xmg(nmg+2:nbg+2,t) = xmg(nmg+1,t) + pmgAcc_s(nmg);
            ymg(nmg+2:nbg+2,t) = ymg(nmg+1,t);
        end
    end

    lbdAcc = importdata(strcat(busL{k},'\tcmpa_lbD.a2m'));
    pbdAcc = importdata(strcat(busL{k},'\tcmpa_pbD.a2m'));
    PDAcc  = importdata(strcat(busL{k},'\tcmpa_PD.a2m'));

    pbdAcc = pbdAcc';
    lbdAcc = lbdAcc';
    PDAcc = PDAcc';

    [nbd,tp] = size(lbdAcc);
    [lbdAcc_s,IXd] = sort(lbdAcc,1,'descend');
    for t = 1:tp
        pbdAcc_s(:,t) = pbdAcc(IXd(:,t),t);
        PDAcc_s(:,t) = PDAcc(IXd(:,t),t);
    end

    xd = zeros(nbd+2,tp); yd = zeros(nbd+2,tp);
    for t = 1:tp
        yd(1,t)=lbdAcc_s(1,t);
        for j = 2:nbd+1
            yd(j,t) = lbdAcc_s(j-1,t);
            if j==2
                xd(j,t) = 0;
            else
                xd(j,t) = xd(j-1,t) + pbdAcc_s(j-2,t);
            end
        end
        xd(nbd+2,t) = xd(nbd+1,t) + pbdAcc_s(nbd,t);
        yd(nbd+2,t) = yd(nbd+1,t);
    end

    xmd = zeros(nbd+2,tp); ymd = zeros(nbd+2,tp);
    for t = 1:tp
        pmdAcc_s = PDAcc_s(PDAcc_s(:,t)>0,t);
        lmdAcc_s = lbdAcc_s(PDAcc_s(:,t)>0,t);
        [nmd,tpm] = size(lmdAcc_s);
        ymd(1,t)=lmdAcc_s(1);
        for j = 2:nmd+1
            ymd(j,t) = lmdAcc_s(j-1);
            if j==2
                xmd(j,t) = 0;
            else
                xmd(j,t) = xmd(j-1,t) + pmdAcc_s(j-2);
            end
        end
        xmd(nmd+2:nbd+2,t) = xmd(nmd+1,t) + pmdAcc_s(nmd);
        ymd(nmd+2:nbd+2,t) = ymd(nmd+1,t);
    end

    MPA_eq = importdata(strcat(busL{k},'\tcmpa_eq.a2m'));
    MPA_eq = MPA_eq';

    NSWT = 0;

    maxenergy = max(xg(nbg+2),xd(nbd+2))+1;
    maxprice1 = max(yg(nbg+2,:));
    maxprice2 = max(yd(1,:));
    maxprice = max(maxprice1,maxprice2)+1;

    epsy = 0.0*maxprice/150;
    ymg = ymg-epsy;
    ymd = ymd+epsy;

    ncol = 2;
    for t=tpS
        energy = MPA_eq(1,t);
        price = MPA_eq(2,t);
        NSW = MPA_eq(3,t);
        NSWT = NSWT + NSW;
        nexttile;
%        subplot(nB,ntp,(k-1)*ntp+t-(tpi-1))
%        subplot(ceil(tp/ncol),ncol,t)
        hold on;
        box on;
        fs = 10;
        if k==min(kS)
            title({['\fontname{arial}{\bf Time period }',num2str(t)];['{\it\lambda}^*=',num2str(price),', {\itP }^*=',num2str(energy),', {\itNSW ^*}= ',num2str(NSW)]},'FontName','Cambria Math','Fontsize',fs,'FontWeight','bold');
        else
            title({['{\it\lambda}^*=',num2str(price),', {\itP }^*=',num2str(energy),', {\itNSW ^*}= ',num2str(NSW)]},'FontName','Arial','Fontsize',fs,'FontWeight','bold');
        end
        if t==min(tpS)
            if nB == 1 
                ylabel('{\it\lambda} [€/MWh]','Fontname','Arial','Fontsize',fs,'FontWeight','bold')
            else
                ylabel({busL{k};'{\it\lambda} [€/MWh]'},'Fontname','Arial','Fontsize',fs,'FontWeight','bold')
            end
        else
%            ylabel('{\it\lambda} [€/MWh]','Fontname','Cambria Math','Fontsize',fs)
        end
        if k==max(kS)
            xlabel('{\itP} [MW]','Fontname','Arial','Fontsize',fs,'FontWeight','bold')
        end
        set(gcf,'Color','w'); set(gca,'Fontsize',fs);

        axis ([0 maxenergy 0 maxprice]);
        stairs(xd(:,t),yd(:,t),'-','Color',[.6,.6,.0],'LineWidth',2)
        stairs(xmd(:,t),ymd(:,t),'-','Color',[0.60,.60,.0],'LineWidth',6)
        stairs(xg(:,t),yg(:,t),'-','Color',[.26,.26,0.8],'LineWidth',2)
        stairs(xmg(:,t),ymg(:,t),'-','Color',[.26,.26,0.8],'LineWidth',6)
        plot([0;energy],[price;price],':k','LineWidth',3);
        plot([energy;energy],[price;0],':k','LineWidth',3);
        plot([energy;energy],[price;price],'xr','MarkerSize',20,'LineWidth',10);
%        hold off;
    end
end



