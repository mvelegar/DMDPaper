% Plot relative reconstruction errors vs Modes in OPTDMD for 40 days of 
% preprocessed data.
%  
% REQUIRES:
%   START and TEND relative error files from results folder

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

R=1:50;
relErrStart=load('../results/relErrVsRStartAllLat30.mat');
relErrStart=relErrStart.relErrStart; 
relErrTend=load('../results/relErrVsRTendAllLat30.mat');
relErrTend=relErrTend.relErrTend; 


%% Compute the relative error and plot
% x is a vector, matrix, or any numeric array of data. NaNs are ignored.
% p is the confidence level (ie, 95 for 95% CI)
% The output is 1x2 vector showing the [lower,upper] interval values.
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));

CIStartAll = NaN(2,length(R),nChems); CITendAll = NaN(2,length(R),nChems);

for r=1:length(R)
    for iChem=1:nChems
        CIStartAll(:,r,iChem)=CIFcn(relErrStart(r,iChem),95);
        CITendAll(:,r,iChem)=CIFcn(relErrTend(r,iChem),95);
    end
end
%% Plot the error bar
fontSize=18;
figure();
ha = tight_subplot(6,2,[.05 .05],[.1 .05],[.1 .05]);

for iChem=1:length(chem_species)
    nPlt=2*iChem-1; axes(ha(nPlt));
    plot(R,relErrStart(:,iChem),'-s','MarkerSize',3,...
                                   'Linewidth',4,'Color',[0 0 0 1]);
    grid on;  xlim([1 length(R)]);
    if mode(nPlt,2) ~= 0 && nPlt ~=11
        set(gca,'LineWidth',4,'FontSize',fontSize,'xtick',...
            5:5:length(R),...
            'xticklabel',[]);
    elseif nPlt==11 
        set(gca,'LineWidth',4,'FontSize',fontSize,'xtick',...
            10:10:length(R),...
            'xticklabel',{'10', '20', '30', '40', '50'});
     end
%--------------------------------------------------------------------------
    nPlt=2*iChem; axes(ha(nPlt));
    plot(R,relErrTend(:,iChem),'-s','MarkerSize',3,...
                                   'Linewidth',4,'Color',[0 0 0 1]);
    grid on;  xlim([1 length(R)]);
    if mod(nPlt,2) == 0 && nPlt ~= 12
        set(gca,'LineWidth',4,'FontSize',fontSize,'xtick',...
            5:5:length(R),...
            'xticklabel',[]);
    elseif nPlt==12
        set(gca,'LineWidth',4,'FontSize',fontSize,'xtick',...
            10:10:length(R),...
            'xticklabel',{'10', '20', '30', '40', '50'});
     end
end
