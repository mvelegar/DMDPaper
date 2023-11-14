% OPTDMD for 2 months of preprocessed data
% Reconstruct and predict the time series using 
% CON2OPTDMD for Lat=30. Save the predicted results.
% Requires: computeSpatialModes results from results folder

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

% The chemical species info
% The 6 chemical species of interest
chem_species=cellstr(...
    ['NO  ';
    'O3  ';
    'NO2 ';
    'CO  ';
    'ISOP';
    'OH  ']);
nChems=length(chem_species);

%% Load the results for OPT & Con2OPT DMD
% Time info
nDaysTotal=61; nSnapsDay = 72; % For snapshots every 20-min
nSnapsTotal=nDaysTotal*nSnapsDay;
nTrainDays=40; nSnapsTrain=nTrainDays*nSnapsDay;
t=linspace(0,nDaysTotal,nSnapsDay*(nDaysTotal)); % Recon+Pred
nPredDays=nDaysTotal-nTrainDays;
nSnapsPred=nSnapsDay*nPredDays;

% Preprocessed data
% size nLon x nSnapsTotal x nChems x nlat
YStartAll = load('../results/YStartLatAll.mat');
YStartAll = YStartAll.YStartLatAll;
YTendAll = load('../results/YTendLatAll.mat');
YTendAll = YTendAll.YTendLatAll;
% Just the training window preprocessed data
% size nLon x nSnapsTrain x nChems x nlat
YStartTrainAll = load('../results/YStartTrainAll.mat');
YStartTrainAll = YStartTrainAll.YStartTrainAll;
YTendTrainAll = load('../results/YTendTrainAll.mat');
YTendTrainAll = YTendTrainAll.YTendTrainAll;

% Reconstructed data from OPT DMD with Con2
% size nLon x nSnapsTrain x nChems x nlat
YStartCon2OPTAll = load('../results/YStartOPTAll.mat');
YStartCon2OPTAll = YStartCon2OPTAll.YStartCon2OPTAll;
YTendCon2OPTAll = load('../results/YTendOPTAll.mat');
YTendCon2OPTAll = YTendCon2OPTAll.YTendCon2OPTAll;
% The constrained OPT DMD results
wStartCon2OPTAll = load('../results/wStartOPTAll.mat');
wStartCon2OPTAll = wStartCon2OPTAll.wStartCon2OPTAll;
bStartCon2OPTAll = load('../results/bStartOPTAll.mat');
bStartCon2OPTAll = bStartCon2OPTAll.bStartCon2OPTAll;
eStartCon2OPTAll = load('../results/eStartOPTAll.mat');
eStartCon2OPTAll = eStartCon2OPTAll.eStartCon2OPTAll;
wTendCon2OPTAll = load('../results/wTendOPTAll.mat');
wTendCon2OPTAll = wTendCon2OPTAll.wTendCon2OPTAll;
eTendCon2OPTAll = load('../results/eTendOPTAll.mat');
eTendCon2OPTAll = eTendCon2OPTAll.eTendCon2OPTAll;
bTendCon2OPTAll = load('../results/bTendOPTAll.mat');
bTendCon2OPTAll = bTendCon2OPTAll.bTendCon2OPTAll;

%% Re-compute the r based on results
% r will chance for iChem and iLat
% Number of reduced modes
rStartOld=[25 25 25 50 25 25];
rTendOld=[20 20 50 20 20 50];

rStartOPT=NaN(nChems,nlat);
rTendOPT=NaN(nChems,nlat);
rStartCon2OPT=NaN(nChems,nlat);
rTendCon2OPT=NaN(nChems,nlat);
for iLat=1:nlat
    for iChem=1:nChems
        if sum(isnan(bStartCon2OPTAll(:,iChem,iLat)))==0
            rStartCon2OPT(iChem,iLat)=length(bStartCon2OPTAll(:,iChem,iLat));
        else
            ind=find(isnan(bStartCon2OPTAll(:,iChem,iLat)),1);
            rStartCon2OPT(iChem,iLat)=ind-1;
        end
        if sum(isnan(bTendCon2OPTAll(:,iChem,iLat)))==0
            rTendCon2OPT(iChem,iLat)=length(bTendCon2OPTAll(:,iChem,iLat));
        else
            ind=find(isnan(bTendCon2OPTAll(:,iChem,iLat)),1);
            rTendCon2OPT(iChem,iLat)=ind-1;
        end
    end
end



%% Part III: Compute prediction results
% sample time
nDaysTest=20; nSnapsTest=nDaysTest*nSnapsDay;
YPrdctStartCon2OPTAll=NaN(length(x),nSnapsTest+1,length(chem_species),nlat);
YPrdctTendCon2OPTAll=NaN(length(x),nSnapsTest+1,length(chem_species),nlat);

for iLat=1:nlat
    for iChem=1:length(chem_species)

        % predicted value based on last training snapshot
        % Here I do not need to compute a new b since this snapshot was already
        % used to compute b in the OPTDMD algorithm

        YPrdctStartCon2OPTAll(:,:,iChem,iLat) = ...
            real(wStartCon2OPTAll(:,1:rStartCon2OPT(iChem,iLat),iChem,iLat)*...
            diag(bStartCon2OPTAll(1:rStartCon2OPT(iChem,iLat),iChem,iLat))*...
            exp(eStartCon2OPTAll(1:rStartCon2OPT(iChem,iLat),iChem,iLat)...
            *t(nSnapsTrain:nSnapsTrain+nSnapsTest)));
        YPrdctTendCon2OPTAll(:,:,iChem,iLat) = ...
            real(wTendCon2OPTAll(:,1:rTendCon2OPT(iChem,iLat),iChem,iLat)*...
            diag(bTendCon2OPTAll(1:rTendCon2OPT(iChem,iLat),iChem,iLat))*...
            exp(eTendCon2OPTAll(1:rTendCon2OPT(iChem,iLat),iChem,iLat)...
            *t(nSnapsTrain:nSnapsTrain+nSnapsTest)));

        % IF I want to use a future snapshot to predict, then use the following
        %     snapFuture = nSnapsTrain;% Last snapshot for training
        %     A_prdct = w*diag(exp(e*Iord(nSnaps+snapFuture)));
        %     b_prdct = A_prdct\YStart(:,nSnaps+snapFuture);
        %     Y3 = w*diag(b_prdct)*exp(e*Iord(nSnaps+1:end)); Y3=real(Y3);



    end
end
% SAVE Results for Error bar plots
save('../results/YPrdctStartCon2OPTAll.mat','YPrdctStartCon2OPTAll'...
    ,'-v7.3');
save('../results/YPrdctTendCon2OPTAll.mat','YPrdctTendCon2OPTAll'...
    ,'-v7.3');

%% Plot the predicted results
lons=-180:5:-155;%:60:150;  %pltLonInds=NaN(size(lons));
iLat=12;
tPlt=t(1:nSnapsTrain+nSnapsTest);
YStartAllPlt=YStartAll(:,1:nSnapsTrain+nSnapsTest,:,iLat);
YTendAllPlt=YTendAll(:,1:nSnapsTrain+nSnapsTest,:,iLat);
fontSize=18;

for iLat=nlat%1:nlat
    for iChem=[1 4]%1:nChems

        % START Results
        % Con2 OPT  Results
        figure();
        ha = tight_subplot(6,1,[.05 .03],[.1 .1],[.1 .05]);
        ymin1=NaN(length(lons),1); ymin2=NaN(length(lons),1);
        ymin3=NaN(length(lons),1);
        ymax1=NaN(length(lons),1); ymax2=NaN(length(lons),1);
        ymax3=NaN(length(lons),1);
        for iLon=1:length(lons)
            ymin1(iLon)=min(YPrdctStartCon2OPTAll(iLon,:,iChem,iLat));
            ymax1(iLon)=max(YPrdctStartCon2OPTAll(iLon,:,iChem,iLat));
            ymin2(iLon)=min(YStartAllPlt(iLon,:,iChem));
            ymax2(iLon)=max(YStartAllPlt(iLon,:,iChem));
            ymin3(iLon)=min(YStartCon2OPTAll(iLon,:,iChem,iLat));
            ymax3(iLon)=max(YStartCon2OPTAll(iLon,:,iChem,iLat));
        end
        yminlim=min(ymin2);%min(min([ymin1 ymin2 ymin3]));
        ymaxlim=max(ymax2);%max(max([ymax1 ymax2 ymax3]));
        yminlim=yminlim-0.001*yminlim; ymaxlim=ymaxlim+0.5*ymaxlim;
        legendStr=cellstr(...
            ['Reconstructed   ';
            'Predicted       ';
            'STARTData       ';
            'PredictionWindow']);
        for iLon= 1:length(lons)
            axes(ha(iLon));
            if iLon==2
                rcnsrctLine = plot(tPlt(1:nSnapsTrain),YStartCon2OPTAll(iLon,:,iChem,iLat),...
                    'm','Linewidth',4); hold on;
                prdctLine = plot(tPlt(nSnapsTrain:end),YPrdctStartCon2OPTAll(iLon,:,iChem,iLat),...
                    'color',[0.8500 0.3250 0.0980],'Linewidth',4); hold on;
                dataLine = plot(tPlt,YStartAllPlt(iLon,:,iChem),':','color',...
                    [0 0 0 0.6],...
                    'Linewidth',4);
                prdWndwLine = plot([tPlt(nSnapsTrain) tPlt(nSnapsTrain)],...
                    [yminlim ymaxlim],'color',[0.4660 0.6740 0.1880],...
                    'Linewidth',4);hold off;
            else
                plot(tPlt(1:nSnapsTrain),YStartCon2OPTAll(iLon,:,iChem,iLat),...
                    'm','Linewidth',4); hold on;
                plot(tPlt(nSnapsTrain:end),YPrdctStartCon2OPTAll(iLon,:,iChem,iLat),...
                    'color',[0.8500 0.3250 0.0980],'Linewidth',4); hold on;
                plot(tPlt,YStartAllPlt(iLon,:,iChem),':','color',[0 0 0 0.5],...
                    'Linewidth',4);
                plot([tPlt(nSnapsTrain) tPlt(nSnapsTrain)],[yminlim ymaxlim],...
                    'color',[0.4660 0.6740 0.1880],'Linewidth',4);hold off;
            end
            grid on;
            xlim([tPlt(nSnapsTrain-10*nSnapsDay) tPlt(end)]);
            ylim([yminlim ymaxlim]);
            if iLon==6
                set(gca,'LineWidth',4,'FontSize',fontSize,'xtick',...
                    tPlt(nSnapsTrain-10*nSnapsDay):5:tPlt(end),...
                    'xticklabel',{'30','35','40','45','50','55','60'});%,...
                %          'ytick',2e6:4e6:10e6);
            else
                set(gca,'LineWidth',4,'FontSize',fontSize,'xtick',...
                    tPlt(nSnapsTrain-10*nSnapsDay):5:tPlt(end),...
                    'xticklabel',[]);%,'ytick',2e6:4e6:10e6);
            end
            if iChem==1
                ylim([0 4e8])
            end
        end
        % Construct a Legend with the data from the sub-plots
        hL = legend([rcnsrctLine,prdctLine,dataLine,prdWndwLine],...
            legendStr,'Orientation','Horizontal'); legend boxoff;
        hL.FontSize = 18; hL.FontName = 'Times New Roman';
        hold off;
        drawnow;

        %TEND Results
        % Con2 OPT  Results
        figure();
        ha = tight_subplot(6,1,[.05 .03],[.1 .1],[.1 .05]);
        ymin1=NaN(length(lons),1); ymin2=NaN(length(lons),1);
        ymin3=NaN(length(lons),1);
        ymax1=NaN(length(lons),1); ymax2=NaN(length(lons),1);
        ymax3=NaN(length(lons),1);
        for iLon=1:length(lons)
            ymin1(iLon)=min(YPrdctTendCon2OPTAll(iLon,:,iChem,iLat));
            ymax1(iLon)=max(YPrdctTendCon2OPTAll(iLon,:,iChem,iLat));
            ymin2(iLon)=min(YTendAllPlt(iLon,:,iChem));
            ymax2(iLon)=max(YTendAllPlt(iLon,:,iChem));
            ymin3(iLon)=min(YTendCon2OPTAll(iLon,:,iChem,iLat));
            ymax3(iLon)=max(YTendCon2OPTAll(iLon,:,iChem,iLat));
        end
        yminlim=min(ymin3);%min(min([ymin1 ymin2 ymin3]));
        ymaxlim=max(ymax3);%max(max([ymax1 ymax2 ymax3]));
        yminlim=yminlim-0.001*yminlim; ymaxlim=ymaxlim+0.5*ymaxlim;
        legendStr=cellstr(...
            ['Reconstructed   ';
            'Predicted       ';
            'TENDData        ';
            'PredictionWindow']);
        for iLon= 1:length(lons)
            axes(ha(iLon));
            if iLon==2
                rcnsrctLine = plot(tPlt(1:nSnapsTrain),YTendCon2OPTAll(iLon,:,iChem,iLat),...
                    'm','Linewidth',4); hold on;
                prdctLine = plot(tPlt(nSnapsTrain:end),YPrdctTendCon2OPTAll(iLon,:,iChem,iLat),...
                    'color',[0.8500 0.3250 0.0980],'Linewidth',4); hold on;
                dataLine = plot(tPlt,YTendAllPlt(iLon,:,iChem),':','color',...
                    [0 0 0 0.6],...
                    'Linewidth',4);
                prdWndwLine = plot([tPlt(nSnapsTrain) tPlt(nSnapsTrain)],...
                    [yminlim ymaxlim],'color',[0.4660 0.6740 0.1880],...
                    'Linewidth',4);hold off;
            else
                plot(tPlt(1:nSnapsTrain),YTendCon2OPTAll(iLon,:,iChem,iLat),...
                    'm','Linewidth',4); hold on;
                plot(tPlt(nSnapsTrain:end),YPrdctTendCon2OPTAll(iLon,:,iChem,iLat),...
                    'color',[0.8500 0.3250 0.0980],'Linewidth',4); hold on;
                plot(tPlt,YTendAllPlt(iLon,:,iChem),':','color',[0 0 0 0.5],...
                    'Linewidth',4);
                plot([tPlt(nSnapsTrain) tPlt(nSnapsTrain)],[yminlim ymaxlim],...
                    'color',[0.4660 0.6740 0.1880],'Linewidth',4);hold off;
            end
            grid on;
            xlim([tPlt(nSnapsTrain-10*nSnapsDay) tPlt(end)]);
            if iLon==6
                set(gca,'LineWidth',4,'FontSize',fontSize,'xtick',...
                    tPlt(nSnapsTrain-10*nSnapsDay):5:tPlt(end),...
                    'xticklabel',{'30','35','40','45','50','55','60'});
            else
                set(gca,'LineWidth',4,'FontSize',fontSize,'xtick',...
                    tPlt(nSnapsTrain-10*nSnapsDay):5:tPlt(end),...
                    'xticklabel',[]);
            end

        end
        % Construct a Legend with the data from the sub-plots
        hL = legend([rcnsrctLine,prdctLine,dataLine,prdWndwLine],...
            legendStr,'Orientation','Horizontal'); legend boxoff;
        hL.FontSize = 18; hL.FontName = 'Times New Roman';
        hold off;
        drawnow;
    end
end

