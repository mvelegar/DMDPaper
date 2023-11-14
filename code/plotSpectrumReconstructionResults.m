% Plot Spectrum and Reconstruction results for 4 DMD algorithms for 2 
% months of preprocessed data for OH chemical species: 
%   (i) plot the spectrum of eigen values from the 4 different DMD algorithms 
%   (ii) plot the reconstruction results for 40 days for 4 different DMD
%   algorithms
% Requires: results folder with results from running compareClassicOPTDMD
% script



clear variables; close all; clc

%% Set up grid, time and species variables
% Landmap information
ncfileMap='../data/landmap.nc';
landmap=ncread(ncfileMap,'LANDMAP');
landmap=permute(landmap,[2 1]);
% Read the Spatial grid
x=ncread(ncfileMap,'lon');% Longitude(-180:5:175)
y=ncread(ncfileMap,'lat');% Latitude(-89,-86:4:86,89)
nLon=length(x); nLat=length(y);
% I am picking one elevation lev=1, all latitudes between
lev=1;
% Latitudes limited such that I do not have to cut out too many snap shots
% to keep day lenghts consistent across a latitude. Right now picking such
% that no more than 4 snap shots are cut off.
latLim=[-14 30];
latVecIndLim(1)=find(y==latLim(1)); latVecIndLim(2)=find(y==latLim(2));
nlat=latVecIndLim(2)-latVecIndLim(1)+1;
yLim=y(latVecIndLim(1):latVecIndLim(2)); % 12 latitudes

% Time info
nSnapsDay = 72; % For snapshots every 20-min
nTrainDays=40; nSnapsTrain=nTrainDays*nSnapsDay;
nPlotDays=15; nSnapsPlot=nPlotDays*nSnapsDay;
% Time vector
t=linspace(0,nTrainDays,nSnapsDay*(nTrainDays)); %in days


% The chemical species info
% The 6 chemical species of interest
chem_species=cellstr(...
    ['NO  ';
    'O3  ';
    'NO2 ';
    'OH  ';
    'ISOP';
    'CO  ';]);

%% Load the results for Spectrum
rStart=[25 25 25 25 25 50];
rTend=[20 20 50 50 20 20]; iLat=12;
eigValsStart=NaN(max(rStart),length(chem_species),4);
eigValsTend=NaN(max(rTend),length(chem_species),4);

for iChem=1:length(chem_species)
    % START data
    omega = load(['../results/omegaStart',chem_species{iChem},...
        'lat30.mat']);
    omega=omega.omega; eigValsStart(1:rStart(iChem),iChem,4)=omega;
    e = load(['../results/eStart',chem_species{iChem},...
        'lat30.mat']);
    e = e.e; eigValsStart(1:rStart(iChem),iChem,1)=e;clear omega e;
    e = load(['../results/eStartCon1',chem_species{iChem},...
        'lat30.mat']);
    e = e.e; eigValsStart(1:rStart(iChem),iChem,3)=e;clear e;
    e = load(['../results/eStartCon2',chem_species{iChem},...
        'lat30.mat']);
    e = e.e; eigValsStart(1:rStart(iChem),iChem,2)=e;clear e;
    % TEND data
    omega = load(['../results/omegaTend',chem_species{iChem},...
        'lat30.mat']);
    omega=omega.omega; eigValsTend(1:rTend(iChem),iChem,4)=omega; 
    e = load(['../results/eTend',chem_species{iChem},...
        'lat30.mat']);
    e = e.e; eigValsTend(1:rTend(iChem),iChem,1)=e;clear omega e;
    e = load(['../results/eTendCon1',chem_species{iChem},...
        'lat30.mat']);
    e = e.e; eigValsTend(1:rTend(iChem),iChem,3)=e; clear e; 
    e = load(['../results/eTendCon2',chem_species{iChem},...
        'lat30.mat']);
    e = e.e; eigValsTend(1:rTend(iChem),iChem,2)=e; clear e;
end
%% Plot the spectrum comparing Classic and OPT DMD
fontSize=18;
% iChem=4 is the OH chemical species
for iChem=4%1:length(chem_species)
    figure();
    ha = tight_subplot(4,2,[.04 .07],[.1 0.1],[.05 .01]);
    temp = eigValsStart(:,iChem,:); temp=squeeze(temp);
    eMinX=min(real(temp(:)));
    eMaxX=max(real(temp(:)));
    eMinY=min(imag(temp(:)));
    eMaxY=max(imag(temp(:)));
    
    iEig=1;
    for i=1:2:7 % The classic DMD spectrum for START
     axes(ha(i)); 
     scatter(real(eigValsStart(:,iChem,iEig)),...
                 imag(eigValsStart(:,iChem,iEig)),75,'k','filled');
     xlim([eMinX eMaxX]); ylim([floor(eMinY) ceil(eMaxY)]); 
     iEig = iEig+1;
     grid on;
     if i==7
       set(gca,'LineWidth',2,'FontSize',fontSize);  
     else
         set(gca,'xticklabel',[],'LineWidth',2,'FontSize',fontSize);
     end
    end
    
    temp = eigValsTend(:,iChem,:); temp=squeeze(temp);
    eMinX=min(real(temp(:)));
    eMaxX=max(real(temp(:)));
    eMinY=min(imag(temp(:)));
    eMaxY=max(imag(temp(:)));
    iEig=1;
    for i=2:2:8 % The classic DMD spectrum for START
     axes(ha(i)); 
     scatter(real(eigValsTend(:,iChem,iEig)),...
                 imag(eigValsTend(:,iChem,iEig)),75,'k','filled');
     xlim([eMinX eMaxX]); ylim([floor(eMinY) ceil(eMaxY)]);
     iEig = iEig+1;
     grid on;
     if i==8
       set(gca,'LineWidth',2,'FontSize',fontSize);  
     else
         set(gca,'xticklabel',[],'LineWidth',2,'FontSize',fontSize);
     end
    end

end


%% Load Data in for the carpet pictures
YStartTrainAll = NaN(length(x),nSnapsTrain,length(chem_species));
YTendTrainAll = NaN(length(x),nSnapsTrain,length(chem_species));
Y1StartAll = NaN(length(x),nSnapsTrain,length(chem_species));
Y2StartAll = NaN(length(x),nSnapsTrain,length(chem_species));
Y3StartAll = NaN(length(x),nSnapsTrain,length(chem_species));
Y4StartAll = NaN(length(x),nSnapsTrain,length(chem_species));
Y1TendAll = NaN(length(x),nSnapsTrain,length(chem_species));
Y2TendAll = NaN(length(x),nSnapsTrain,length(chem_species));
Y3TendAll = NaN(length(x),nSnapsTrain,length(chem_species));
Y4TendAll = NaN(length(x),nSnapsTrain,length(chem_species));


for iChem = 1:length(chem_species)
    Y = load(['../results/YStartTrain',chem_species{iChem},'.mat']);
    Y = Y.YStartTrain; YStartTrainAll(:,:,iChem)=real(Y); clear Y;
    Y = load(['../results/Y1Start',chem_species{iChem},'.mat']);
    Y = Y.Y1; Y1StartAll(:,:,iChem)=real(Y); clear Y;
    Y = load(['../results/Y2Start',chem_species{iChem},'.mat']);
    Y = Y.Y2; Y2StartAll(:,:,iChem)=real(Y); clear Y;
    Y = load(['../results/Y3Start',chem_species{iChem},'.mat']);
    Y = Y.Y3; Y3StartAll(:,:,iChem)=real(Y); clear Y;
    Y = load(['../results/Y4Start',chem_species{iChem},'.mat']);
    Y = Y.Y4; Y4StartAll(:,:,iChem)=real(Y); clear Y;
    
    Y = load(['../results/YTendTrain',chem_species{iChem},'.mat']);
    Y = Y.YTendTrain; YTendTrainAll(:,:,iChem)=real(Y);
    Y = load(['../results/Y1Tend',chem_species{iChem},'.mat']);
    Y = Y.Y1; Y1TendAll(:,:,iChem)=real(Y); clear Y;
    Y = load(['../results/Y2Tend',chem_species{iChem},'.mat']);
    Y = Y.Y2; Y2TendAll(:,:,iChem)=real(Y); clear Y;
    Y = load(['../results/Y3Tend',chem_species{iChem},'.mat']);
    Y = Y.Y3; Y3TendAll(:,:,iChem)=real(Y); clear Y;
    Y = load(['../results/Y4Tend',chem_species{iChem},'.mat']);
    Y = Y.Y4; Y4TendAll(:,:,iChem)=real(Y); clear Y;
end


%% Plot the carpet pictures


% For the plots
addpath('cmocean_v1.4/cmocean/');
fontSize=18;

% iChem=4 is the OH chemical species
for iChem=4%1:length(chem_species)
    YStartTrain = YStartTrainAll(:,1:nSnapsPlot,iChem);
    YTendTrain = YTendTrainAll(:,1:nSnapsPlot,iChem);
    Y1Start = Y1StartAll(:,1:nSnapsPlot,iChem);
    Y2Start = Y2StartAll(:,1:nSnapsPlot,iChem);
    Y3Start = Y3StartAll(:,1:nSnapsPlot,iChem);
    Y4Start = Y4StartAll(:,1:nSnapsPlot,iChem);
    Y1Tend = Y1TendAll(:,1:nSnapsPlot,iChem);
    Y2Tend = Y2TendAll(:,1:nSnapsPlot,iChem);
    Y3Tend = Y3TendAll(:,1:nSnapsPlot,iChem);
    Y4Tend = Y4TendAll(:,1:nSnapsPlot,iChem);
    figure();
    
    ha = tight_subplot(5,2,[.01 .07],[.1 0.01],[.05 .1]);
    cmin1=min(YStartTrain(:)); cmin2=min(Y1Start(:)); cmin3=min(Y2Start(:));
    cmin4=min(Y3Start(:)); cmin5=min(Y4Start(:));
    cminStart=min([cmin1 cmin2 cmin3 cmin4 cmin5]);
    cmax1=max(YStartTrain(:)); cmax2=max(Y1Start(:)); cmax3=max(Y2Start(:));
    cmax4=max(Y3Start(:)); cmax5=max(Y4Start(:));
    cmaxStart=max([cmax1 cmax2 cmax3 cmax4 cmax5]);
    cmin1=min(YTendTrain(:)); cmin2=min(Y1Tend(:)); cmin3=min(Y2Tend(:));
    cmin4=min(Y3Tend(:)); cmin5=min(Y4Tend(:));
    cminTend=min([cmin1 cmin2 cmin3 cmin4 cmin5]);
    cmax1=max(YTendTrain(:)); cmax2=max(Y1Tend(:)); cmax3=max(Y2Tend(:));
    cmax4=max(Y3Tend(:)); cmax5=max(Y4Tend(:));
    cmaxTend=max([cmax1 cmax2 cmax3 cmax4 cmax5]);
    
    tPlot=t(1:nSnapsPlot);
    cmap = cmocean('balance');
    axes(ha(1)); % Data
    pcolor(tPlot,x,YStartTrain);
    set(gca,...
        'ytick',[-90,0,90],'LineWidth',2,...
        'FontSize',fontSize,'YTickLabel',...
        ({'-90','0','90'}));
    colormap(cmap);
    shading interp;
    % ylabel('$lon$','Interpreter','Latex','Fontsize',fontSize);
    %     caxis([cminStart cmaxStart]);
    axes(ha(9)); %DMD
    pcolor(tPlot,x,Y1Start);
    set(gca,...
        'ytick',[-90,0,90],'LineWidth',2,...
        'FontSize',fontSize,'YTickLabel',...
        ({'-90','0','90'}));
    colormap(cmap);
    shading interp;
    % ylabel('$lon$','Interpreter','Latex','Fontsize',fontSize);
    caxis([cminStart cmaxStart]);
    axes(ha(3));% OPTDMD
    pcolor(tPlot,x,Y2Start);
    set(gca,'xtick',[0,10,20,30,40],...
        'XTickLabel',...
        ({'0','10','20','30','40'}),'LineWidth',2,...
        'ytick',[-90,0,90],...
        'FontSize',fontSize,'YTickLabel',...
        ({'-90','0','90'}));
    colormap(hsv);
    shading interp;
    
    caxis([cminStart cmaxStart]);
    axes(ha(7));
    pcolor(tPlot,x,Y3Start);
    set(gca,'xtick',[0,10,20,30,40],...
        'XTickLabel',...
        ({'0','10','20','30','40'}),'LineWidth',2,...
        'ytick',[-90,0,90],...
        'FontSize',fontSize,'YTickLabel',...
        ({'-90','0','90'}));
    colormap(cmap);
    shading interp;
    caxis([cminStart cmaxStart]);
    axes(ha(5));
    pcolor(tPlot,x,Y4Start);
    set(gca,'xtick',[0,10,20,30,40],...
        'XTickLabel',...
        ({'0','10','20','30','40'}),'LineWidth',2,...
        'ytick',[-90,0,90],...
        'FontSize',fontSize,'YTickLabel',...
        ({'-90','0','90'}));
    colormap(cmap);
    shading interp;
    caxis([cminStart cmaxStart]);
    hp1 = get(ha(1),'Position'); hp9 = get(ha(9),'Position');
    cb_x = hp9(1) + hp9(3) + 0.005;
    cb_y = hp9(2);
    cb_w = 0.025;
    cb_h = hp1(2) + hp1(4) - hp9(2);
    colorbar('Position', [cb_x cb_y cb_w cb_h]);
    drawnow;
    % ------------Tend ones----------------------------------------------------
    cmap = cmocean('balance');
    axes(ha(2));
    pcolor(tPlot,x,YTendTrain);
    set(gca,...
        'ytick',[],...%[-90,0,90],...
        'LineWidth',2,...
        'FontSize',fontSize);
    % center the colormap
    cLimCenterTend=max(abs([cminTend cmaxTend]));
    caxis([-cLimCenterTend cLimCenterTend]);
    colormap(cmap);
    shading interp;
    % ylabel('$lon$','Interpreter','Latex','Fontsize',fontSize);
    
    axes(ha(10));
    pcolor(tPlot,x,Y1Tend);
    set(gca,...
        'ytick',[],...%[-90,0,90],
        'LineWidth',2,'FontSize',fontSize);
    % center the colormap
    caxis([-cLimCenterTend cLimCenterTend]);
    colormap(cmap);
    shading interp;
    % ylabel('$lon$','Interpreter','Latex','Fontsize',fontSize);
    
    axes(ha(4));
    pcolor(tPlot,x,Y2Tend);
    set(gca,'xtick',[0,10,20,30,40],...
        'LineWidth',2,...
        'ytick',[],...%[-90,0,90],...
        'FontSize',fontSize);
    % center the colormap
    caxis([-cLimCenterTend cLimCenterTend]);
    colormap(cmap);
    shading interp;
    axes(ha(8));
    pcolor(tPlot,x,Y3Tend);
    set(gca,'xtick',[0,10,20,30,40],...
        'LineWidth',2,...
        'ytick',[],...%[-90,0,90],...
        'FontSize',fontSize);
    % center the colormap
    caxis([-cLimCenterTend cLimCenterTend]);
    colormap(cmap);
    shading interp;
    axes(ha(6));
    pcolor(tPlot,x,Y4Tend);
    set(gca,'xtick',[0,10,20,30,40],...
        'LineWidth',2,...
        'ytick',[],...%[-90,0,90],...
        'FontSize',fontSize);
    % center the colormap
    caxis([-cLimCenterTend cLimCenterTend]);
    colormap(cmap);
    shading interp;
    set(ha(1:10),'XTickLabel','','XTick',[]);
    hp2 = get(ha(2),'Position'); hp10 = get(ha(10),'Position');
    cb_x = hp10(1) + hp10(3) + 0.005;
    cb_y = hp10(2);
    cb_w = 0.025;
    cb_h = hp2(2) + hp2(4) - hp10(2);
    colorbar('Position', [cb_x cb_y cb_w cb_h]);
    
end