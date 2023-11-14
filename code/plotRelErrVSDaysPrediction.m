% Section 4
% Compute and plot relative error in prediction for 20 days using 
% CON2OPTDMD for lat 30
% Requires: computeSpatialModes and computePlotTimeSeriesReconstructPredict 
% results from results folder

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
    'OH  ';
    'ISOP';
    'CO  ';]);
nChems=length(chem_species);

%% Load the results for OPT & Con2OPT DMD
% Time info
nDaysTotal=61; nSnapsDay = 72; % For snapshots every 20-min
nSnapsTotal=nDaysTotal*nSnapsDay;
nTrainDays=40; nSnapsTrain=nTrainDays*nSnapsDay;
t=linspace(0,nDaysTotal,nSnapsDay*(nDaysTotal)); % Recon+Pred
nPredDays=nDaysTotal-nTrainDays;
nSnapsPred=nSnapsDay*nPredDays;
nDaysTest=20; nSnapsTest=nDaysTest*nSnapsDay;

% Preprocessed data
% size nLon x nSnapsTotal x nChems x nlat
YStartAll = load('../results/YStartLatAll.mat');
YStartAll = YStartAll.YStartLatAll;
YTendAll = load('../results/YTendLatAll.mat');
YTendAll = YTendAll.YTendLatAll;

% The prediction results for 21 days
% size nlon x nTestDays(20)*nSnapsDay+1 x nChems x nlat
YPrdctStartCon2OPTAll=load('../results/YPrdctStartCon2OPTAll.mat');
YPrdctStartCon2OPTAll=YPrdctStartCon2OPTAll.YPrdctStartCon2OPTAll;
YPrdctTendCon2OPTAll=load('../results/YPrdctTendCon2OPTAll.mat');
YPrdctTendCon2OPTAll=YPrdctTendCon2OPTAll.YPrdctTendCon2OPTAll;

%% Compute the relative error and plot
% x is a vector, matrix, or any numeric array of data. NaNs are ignored.
% p is the confidence level (ie, 95 for 95% CI)
% The output is 1x2 vector showing the [lower,upper] interval values.
%CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
CIFcn = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),sum(~isnan(x(:)))-1) + mean(x(:),'omitnan'); 

YStartTestAll = YStartAll(:,nSnapsTrain+1:nSnapsTrain+nSnapsTest,:,:);
YTendTestAll = YTendAll(:,nSnapsTrain+1:nSnapsTrain+nSnapsTest,:,:);

CIDaysStartCon2OPTAll = NaN(2,nDaysTest,length(chem_species),nlat);
relErrDayStartCon2OPTAll = NaN(nDaysTest,length(chem_species),nlat);

CIDaysTendCon2OPTAll = NaN(2,nDaysTest,length(chem_species),nlat);
relErrDayTendCon2OPTAll = NaN(nDaysTest,length(chem_species),nlat);

for iLat = 1:nlat
for iChem = 1:length(chem_species)
    sclVleStart = norm(YStartTestAll(:,:,iChem,iLat));
    sclVleTend = norm(YTendTestAll(:,:,iChem,iLat));
    for iDay = 1:nDaysTest
        relErr=(YStartTestAll(:,((iDay-1)*nSnapsDay+1:...
                                         (iDay)*nSnapsDay),iChem,iLat)-...
                        YPrdctStartCon2OPTAll(:,((iDay-1)*nSnapsDay+1:...
                                       (iDay)*nSnapsDay),iChem,iLat))/...
                            sclVleStart;
        relErrDay=mean(relErr); %average it over all 72 lons
        CIDaysStartCon2OPTAll(:,iDay,iChem,iLat)=CIFcn(relErrDay,95); %CI = CIFcn(x,95);
        relErrDayStartCon2OPTAll(iDay,iChem,iLat)=mean(relErrDay); %average over all snap shots in a day

        relErr=(YTendTestAll(:,((iDay-1)*nSnapsDay+1:...
                                         (iDay)*nSnapsDay),iChem,iLat)-...
                        YPrdctTendCon2OPTAll(:,((iDay-1)*nSnapsDay+1:...
                                       (iDay)*nSnapsDay),iChem,iLat))/...
                                       sclVleTend;
        relErrDay=mean(relErr); %average it over all 72 lons
        CIDaysTendCon2OPTAll(:,iDay,iChem,iLat)=CIFcn(relErrDay,95); %CI = CIFcn(x,95);
        relErrDayTendCon2OPTAll(iDay,iChem,iLat)=mean(relErrDay); %average over all snap shots in a day
    end
end
end

%% Plot the error bar
fontSize=18;
figure();
ha = tight_subplot(6,2,[.05 .05],[.1 .05],[.1 .05]);
for iChem=1:length(chem_species)
    nPlt=2*iChem-1; axes(ha(nPlt));
    errorbar(1:nDaysTest,relErrDayStartCon2OPTAll(:,iChem,iLat),CIDaysStartCon2OPTAll(1,:,iChem,iLat),...
             CIDaysStartCon2OPTAll(2,:,iChem,iLat),'-s','MarkerSize',5,...
                                   'MarkerEdgeColor','red',...
                                   'MarkerFaceColor',[1 .6 .6],...
                                   'Linewidth',4,'Color',[0 0 0 0.5]);
    grid on;  xlim([1 nDaysTest]);
    if nPlt==11
         set(gca,'LineWidth',4,'FontSize',fontSize,'xtick',...
            5:5:nDaysTest,...
            'xticklabel',{'5','10','15','20'});
    else
        set(gca,'LineWidth',4,'FontSize',fontSize,'xtick',...
            5:5:nDaysTest,...
            'xticklabel',[]);
    end
    nPlt=2*iChem; axes(ha(nPlt));
    errorbar(1:nDaysTest,relErrDayTendCon2OPTAll(:,iChem,iLat),CIDaysTendCon2OPTAll(1,:,iChem,iLat),...
             CIDaysTendCon2OPTAll(2,:,iChem,iLat),'-s','MarkerSize',5,...
                                   'MarkerEdgeColor','red',...
                                   'MarkerFaceColor',[1 .6 .6],...
                                   'Linewidth',4,'Color',[0 0 0 0.5]);
    grid on;  xlim([1 nDaysTest]);
     if mod(nPlt,2) == 0 && nPlt ~= 12
        set(gca,'LineWidth',4,'FontSize',fontSize,'xtick',...
            5:5:nDaysTest,...
            'xticklabel',[]);
    elseif nPlt==11 || nPlt==12
        set(gca,'LineWidth',4,'FontSize',fontSize,'xtick',...
            5:5:nDaysTest,...
            'xticklabel',{'5','10','15','20'});
    else
        set(gca,'LineWidth',4,'FontSize',fontSize,'xtick',...
            5:5:nDaysTest,...
            'xticklabel',[]);
     end
end

