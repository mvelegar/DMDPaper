% Generate the spatial modes for all latitudes and save them for plotting
% Requires: Original matrices from data folder
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
% Time info
nDaysTotal=61; nSnapsDay = 72; % For snapshots every 20-min
nSnapsTotal=nDaysTotal*nSnapsDay;
nTrainDays=40; nSnapsTrain=nTrainDays*nSnapsDay;
t=linspace(0,nDaysTotal,nSnapsDay*(nDaysTotal)); % Recon+Pred
nPredDays=nDaysTotal-nTrainDays;
nSnapsPred=nSnapsDay*nPredDays;
nPlotDays=15; nSnapsPlot=nPlotDays*nSnapsDay;
% Time vector
tTrain=linspace(0,nTrainDays,nSnapsDay*(nTrainDays)); %in days



% The chemical species info
% The 6 chemical species of interest
chem_species=cellstr(...
    ['NO  ';
    'O3  ';
    'NO2 ';
    'OH  ';
    'ISOP';
    'CO  ';]);
%% Isolate the land cells from ocean cells
landMapLim = landmap(latVecIndLim(1):latVecIndLim(2),:);


%% Load the preprocessed data
% Shifted Raw data 
AStartShift3Months = load('../data/AStartShift3Months.mat');
AStartShift3Months = AStartShift3Months.AStartShift3Months;

ATendShift3Months = load('../data/ATendShift3Months.mat');
ATendShift3Months = ATendShift3Months.ATendShift3Months;

%% Load the results for Spectrum
rStart=[25 25 25 25 25 50];
rTend=[20 20 50 50 20 20]; 
%% Part I: comparison of DMD and OPTDMD
% Fix number of modes r=30. We are only doing a reconstruction for
% fixed lev=1, fixed lat=30, and fixed chem=O3. This is the pretty carpet
% picture
clc; close all;
nTrainDays=40; nPredDays=5; 
rStart=[25 25 25 25 25 50]; 
rTend=[20 20 50 50 20 20]; 

% Observations:
% For iChem = 3; r=50 does a better job for Start
% For iChem = 4; r=40 does a better job for Tend
YStartLatAll=NaN(nLon,nSnapsTotal,length(chem_species),nlat);
YTendLatAll=NaN(nLon,nSnapsTotal,length(chem_species),nlat);

YStartTrainLatAll=NaN(nLon,nSnapsTrain,length(chem_species),nlat);
YTendTrainLatAll=NaN(nLon,nSnapsTrain,length(chem_species),nlat);

YStartCon2OPTAll = NaN(nLon,nSnapsTrain,length(chem_species),nlat);
YTendCon2OPTAll = NaN(nLon,nSnapsTrain,length(chem_species),nlat);
wStartCon2OPTAll = NaN(nLon,max(rStart),length(chem_species),nlat);
wTendCon2OPTAll = NaN(nLon,max(rTend),length(chem_species),nlat);
bStartCon2OPTAll = NaN(max(rStart),length(chem_species),nlat);
bTendCon2OPTAll = NaN(max(rTend),length(chem_species),nlat);
eStartCon2OPTAll = NaN(max(rStart),length(chem_species),nlat);
eTendCon2OPTAll = NaN(max(rTend),length(chem_species),nlat);
for iLat=1:nlat
    for iChem=1:length(chem_species)
    % Original data: time x nLon x nlat x chemSpecies
    YStart=AStartShift3Months(:,:,iLat,iChem);
    YTend=ATendShift3Months(:,:,iLat,iChem);
    
    % Now Dim 1: Lon, Dim 2: time
    YStart=permute(YStart,[2 1]);
    YTend=permute(YTend,[2 1]);
    
    % Take only the training data
    YStartTrain=YStart(:,1:nSnapsTrain);
    YTendTrain=YTend(:,1:nSnapsTrain);
    
    YStartLatAll(:,:,iChem,iLat)=YStart(:,1:nSnapsTotal);
    YTendLatAll(:,:,iChem,iLat)=YTend(:,1:nSnapsTotal);
    YStartTrainLatAll(:,:,iChem,iLat)=YStartTrain;
    YTendTrainLatAll(:,:,iChem,iLat)=YTendTrain;
    
    
    %% OPTDMD with Linear Constraints
    % Start
    lbc = [-Inf*ones(rStart(iChem),1); -Inf*ones(rStart(iChem),1)];
    ubc = [zeros(rStart(iChem),1); Inf*ones(rStart(iChem),1)];
    imode = 1;
    copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);
    opts=varpro_opts('maxiter',300);
    [w,e,b] = optdmd(YStartTrain,t(1:nSnapsTrain),rStart(iChem),imode,opts,[],[],copts);
    clc;
    wStartCon2OPTAll(:,1:rStart(iChem),iChem,iLat)=w;
    eStartCon2OPTAll(1:rStart(iChem),iChem,iLat)=e;
    bStartCon2OPTAll(1:rStart(iChem),iChem,iLat)=b;
    
    % reconstructed values
    Y2 = w*diag(b)*exp(e*t(1:nSnapsTrain)); Y2Start=real(Y2);
    YStartCon2OPTAll(:,:,iChem,iLat)=Y2Start;
    
    % TEND 
    lbc = [-Inf*ones(rTend(iChem),1); -Inf*ones(rTend(iChem),1)];
    ubc = [zeros(rTend(iChem),1); Inf*ones(rTend(iChem),1)];
     copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);
    [w,e,b] = optdmd(YTendTrain,t(1:nSnapsTrain),rTend(iChem),imode,...
              opts);
    clc;
    wTendCon2OPTAll(:,1:rTend(iChem),iChem,iLat)=w;
    eTendCon2OPTAll(1:rTend(iChem),iChem,iLat)=e;
    bTendCon2OPTAll(1:rTend(iChem),iChem,iLat)=b;
    % reconstructed values
    Y2 = w*diag(b)*exp(e*t(1:nSnapsTrain)); Y2Tend=real(Y2);
    YTendCon2OPTAll(:,:,iChem,iLat)=Y2Tend;
    
    end
end
save('../results/YStartLatAll.mat','YStartLatAll','-v7.3');
save('../results/YTendLatAll.mat','YTendLatAll','-v7.3');
save('../results/YStartTrainAll.mat','YStartTrainAll','-v7.3');
save('../results/YTendTrainAll.mat','YTendTrainAll','-v7.3');
save('../results/wStartOPTAll.mat','wStartCon2OPTAll','-v7.3');
save('../results/eStartOPTAll.mat','eStartCon2OPTAll','-v7.3');
save('../results/bStartOPTAll.mat','bStartCon2OPTAll','-v7.3');
save('../results/YStartOPTAll.mat','YStartCon2OPTAll','-v7.3');
save('../results/YTendOPTAll.mat','YTendCon2OPTAll','-v7.3'); 
save('../results/wTendOPTAll.mat','wTendCon2OPTAll','-v7.3');
save('../results/eTendOPTAll.mat','eTendCon2OPTAll','-v7.3');
save('../results/bTendOPTAll.mat','bTendCon2OPTAll','-v7.3');

