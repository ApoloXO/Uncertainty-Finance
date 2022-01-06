 %% Main file: UT_sharedCode.m
% Uncertain Finance
% European options risks by uncertainty theory
% created on Jan 1 2022 
% @author: Carlos Alexander Grajales
% Universidad de Antioquia, Medellín, Colombia
% alexander.grajales@udea.edu.co
% **********************************************************
%
% See working paper:  
% Uncertainty and stochastic theories on European options valuation
% ... and their delta and vega risks.
% Grajales - Medina - Mongrut, January 1, 2022 
% Available at SSRN: https://ssrn.com/abstract=3998282
% **********************************************************

% Stochastic model, BSM: dX_t = r X_t dt + sigma X_t dW_t (1)
% Uncertain model, Liu: dX_t = e X_t dt + sigma X_t dC_t  (2)
% initial value X(0)=X_0
%
% (a) Estimations of stock options: [sc,sp], [uc,up] 
%     [sc,sp] stochastic call, stochastic put - from BSM model
%     [uc,up] uncertain call, uncertain put - from Liu model
% (b) Comparisons of stock options
%     -norm metrics on matrices [uc_sc, up_sp]: Frobenius, norm1, norminf
% (c) follow (a) and (b) for risk letters delta and vega
% 
% for paper we set:
% we examine (a) and (b) at t = 1
% X_0=1; 
% grid sizes for e, sigma, K, T: nume=10; numsigma=10; numk=7; numt=4;


clear;clc;

%% Inputs
X_0=1;
t=1;
Rt=0.01;

nume=10;
numsigma=10;
v_e=linspace(0.01,0.30,nume); 
v_sigma=linspace(0.10,0.60,numsigma);
[sigma,e]=meshgrid(v_sigma,v_e);

numk=7;
numt=4;
Str=linspace(0.7*X_0,1.3*X_0,numk); 
Tm=linspace(0.25*t,1*t,numt);

%% stock options prices [sc,sp], [uc,up] 

% stochastic stock options
% we take variables: Strike, Time, Volatility (K,T,VOL)

Vol=v_sigma; 
[sStr,sTm,sVol]=ndgrid(Str,Tm,Vol); %3D arrange as (K x T) x VOL
[sc,sp]=blsprice(X_0,sStr,Rt,sTm,sVol); % written as (7x40) matrix
sc=reshape(sc,numk,numt,numsigma);  % written as (7x4x10) matrix
sp=reshape(sp,numk,numt,numsigma);

% uncertain stocks options
% we take variables: Strike, Time, Mean, Volatility (K,T,Mu,VOL)

Mn=v_e;
[uStr,uTm,uMn,uVol]=ndgrid(Str,Tm,Mn,Vol); %4D ndgrid
Jc=@(x,e,sigma,Y_0,K,s)max((Y_0*exp(e.*s+(sigma.*s.*sqrt(3)./pi).*(log(x)-log(1-x))))-K,0);
uc=exp(-Rt.*uTm).*integral(@(x)Jc(x,uMn,uVol,X_0,uStr,uTm),0,1,'ArrayValued',true);
Jp=@(x,e,sigma,Y_0,K,s)max(K-(Y_0*exp(e.*s+(sigma.*s.*sqrt(3)./pi).*(log(x)-log(1-x)))),0);
up=exp(-Rt.*uTm).*integral(@(x)Jp(x,uMn,uVol,X_0,uStr,uTm),0,1,'ArrayValued',true);

%% Figures

% call BSM vs LIU, for each vidx, see effects in all the eidx 
for vidx=1:10
    figure(1); clf
    hold on;
    h2=surf(sStr(:,:,vidx),sTm(:,:,vidx),sc(:,:,vidx),'FaceAlpha',0.7); % example of K-T for v_sigma=v_sigma(vidx)
    for eidx=1:10
        h1=mesh(uStr(:,:,eidx,vidx),uTm(:,:,eidx,vidx),uc(:,:,eidx,vidx),'FaceAlpha',0.3); % example of K-T for v_sigma=v_sigma(vidx), v_e=v_e(eidx)
    end
    xlabel('K'); ylabel('T');zlabel('call price')
    title('call price');
    legend([h1, h2], {'uncertain (uc)', 'stochastic (sc)'},'Location','northeast');
    axis tight;
    colormap(winter)
    view(-30,40)
    grid on
    box on
    hold off
    pause
end

% put BSM vs LIU, for each vidx, see effects in all the eidx
for vidx=1:10
    figure(2); clf
    hold on;
    h1=surf(sStr(:,:,vidx),sTm(:,:,vidx),sp(:,:,vidx),'FaceAlpha',0.5); % example of K-T for v_sigma=v_sigma(vidx)
    for eidx=1:10
        h2=mesh(uStr(:,:,eidx,vidx),uTm(:,:,eidx,vidx),up(:,:,eidx,vidx),'FaceAlpha',0.3); % example of K-T for v_sigma=v_sigma(vidx), v_e=v_e(eidx)
    end
    xlabel('K'); ylabel('T');zlabel('put price')
    title('put price');
    legend([h1, h2], {'stochastic (sp)', 'uncertain (up)'},'Location','northeast');
    axis tight;
    colormap(winter)
    view(-150,30)
    grid on
    box on
    hold off
    pause
end

