%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOP-DMD implementation for atmospheric data 
% I have 61 days of data, 40 days of training and 20 days of testing
% Requires: Train files from results folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% The 2 chemical species of interest
% Pick NO and CO like in predicted time series results
chem_species=cellstr(...
    ['NO  ';
    'OH  ';]);
nChems=length(chem_species);

%% Load the data, set up time info
% Time info
nDaysTotal=89; nSnapsDay = 72; % For snapshots every 20-min
nSnapsTotal=nDaysTotal*nSnapsDay;
% Pick 2/3 snapshots for training(60 days) and the rest for testing(29
% days)
nTrainDays=60; nSnapsTrain=nTrainDays*nSnapsDay;
t=linspace(0,nDaysTotal,nSnapsDay*(nDaysTotal)); % Recon+Pred
tTrain=t(1:nSnapsTrain); tTest=t(nSnapsTrain+1:end);
nTestDays=nDaysTotal-nTrainDays;
nSnapsTest=nSnapsDay*nTestDays;

% Load Data in
% Preprocessed data
YStartAll =  NaN(length(x),nSnapsTotal,length(chem_species));
YTendAll =  NaN(length(x),nSnapsTotal,length(chem_species));

for iChem = 1:length(chem_species)
    % Data for the 3 months (July+Aug+Sep) 89 days)
    Y = load(['../results/YStart',chem_species{iChem},'.mat']);
    Y = Y.YStart; YStartAll(:,:,iChem)=real(Y); clear Y;
    Y = load(['../results/YTend',chem_species{iChem},'.mat']);
    Y = Y.YTend; YTendAll(:,:,iChem)=real(Y); clear Y;
end
% Seperate into Train and Test data
YStartAllTrain=YStartAll(:,1:nSnapsTrain,:);
YTendAllTrain=YTendAll(:,1:nSnapsTrain,:);

YStartAllTest=YStartAll(:,nSnapsTrain+1:end,:);
YTestAllTest=YTendAll(:,nSnapsTrain+1:end,:);

%% The BOP-DMD, for 2 chemicals NO and OH

% Number of modes
numModes=9; % discuss this with Nathan. 5 unique b-eval-evec pairs.

% Batch Size p. I have a total of 60 days(nSnapsTrain=4320snapshots). 
% I am using a batch size of 3 days (5%) of my total snapshots(144). 
p = 3*nSnapsDay;

% Number of cycles K, using 100 cycles for now
numCycles = 100; 

% Set up the metrices for lambda/b/w
lambdaStartEnsembleDMD=NaN(numModes,numCycles,length(chem_species));
bStartEnsembleDMD=NaN(numModes,numCycles,length(chem_species));
wStartEnsembleDMD=NaN(nLon,numModes*numCycles,length(chem_species));
lambdaTendEnsembleDMD=NaN(numModes,numCycles,length(chem_species));
bTendEnsembleDMD=NaN(numModes,numCycles,length(chem_species));
wTendEnsembleDMD=NaN(nLon,numModes*numCycles,length(chem_species));


% Set up constraints for OPTDMD (eVals in left-half plane)
lbc = [-Inf*ones([numModes,1]); -Inf*ones([numModes,1])];
ubc = [zeros([numModes,1]); Inf*ones([numModes,1])];
copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);
opts=varpro_opts('maxiter',300,'ifprint',0);

for iChem=1:length(chem_species) 
    % Initialize BOP-DMD; compute e0/b0/w0.
    [wOptStart, eOptStart, bOptStart] = optdmd(YStartAllTrain(:,:,iChem),...
                                       tTrain,numModes,2,opts,[],[],copts);
    
    [wOptTend, eOptTend, bOptTend] = optdmd(YTendAllTrain(:,:,iChem),...
                                       tTrain,numModes,2,opts,[],[],copts);                               
    % The BOP-DMD cycles
    for iCycle=1:numCycles
        %randomly  select p indices from nSnapsTrain
        pInds = randperm(nSnapsTrain,p);
        %sort these p indices
        pInds = sort(pInds);
        
        %Create time and data vectors for this cycle
        tCycle = tTrain(pInds);
        YStartCycle = YStartAllTrain(:,pInds,iChem);
        YTendCycle = YTendAllTrain(:,pInds,iChem);
        
        %START OPTDMD for this cycle with eOPT as initial conditions
        [wStartCycle,eStartCycle,bStartCycle] = optdmd(YStartCycle,tCycle,...
                                    numModes,2,opts,eOptStart,[],copts);
        %sort by imaginary component of evals
        [~,sortImagInd] = sort(imag(eStartCycle));
        
        %save the model generated for each cycle
        wStartEnsembleDMD(:,(iCycle-1)*numModes+1:iCycle*numModes,iChem)=...
                                                wStartCycle(:,sortImagInd);
        lambdaStartEnsembleDMD(:,iCycle,iChem) = eStartCycle(sortImagInd);
        bStartEnsembleDMD(:,iCycle,iChem)=bStartCycle(sortImagInd);
        
        %TEND OPTDMD for this cycle with eOPT as initial conditions
        [wTendCycle,eTendCycle,bTendCycle] = optdmd(YTendCycle,tCycle,...
                                    numModes,2,opts,eOptTend,[],copts);
        %sort by imaginary component of evals
        [~,sortImagInd] = sort(imag(eTendCycle));
        
        %save the model generated for each cycle
        wTendEnsembleDMD(:,(iCycle-1)*numModes+1:iCycle*numModes,iChem)=...
                                                wTendCycle(:,sortImagInd);
        lambdaTendEnsembleDMD(:,iCycle,iChem) = eTendCycle(sortImagInd);
        bTendEnsembleDMD(:,iCycle,iChem)=bTendCycle(sortImagInd);
        clc;
    end
                                               
end


%% Temporal uncertainty quantification: plot the histogram of eigen values

close all; clc; fontSize=18;
for iChem=1:length(chem_species)
    figure(); ha = tight_subplot(round(numModes/2),2,[.05 .1],[.1 0.05],[.05 .07]);   
    for iMode=1:round(numModes/2)
        axes(ha(2*(iMode-1)+1));
        histfit(abs(lambdaStartEnsembleDMD(iMode,:,iChem)),50);
        set(gca,'LineWidth',4,'FontSize',fontSize)
        
        axes(ha(2*(iMode)));
        histfit(abs(lambdaTendEnsembleDMD(iMode,:,iChem)),50);
        set(gca,'LineWidth',4,'FontSize',fontSize)
    end  
    
end

clc;
% Remove outliers for trimmed data - 10% threshold
for iChem=1:length(chem_species)
figure(); ha = tight_subplot(round(numModes/2),2,[.05 .1],[.1 0.05],[.05 .07]);    
    for iMode=1:round(numModes/2)
        axes(ha(2*(iMode-1)+1));
        histfit(rmoutliers(abs(lambdaStartEnsembleDMD(iMode,:,iChem)),...
            "percentiles",[10 90]),50);
        set(gca,'LineWidth',4,'FontSize',fontSize)

        axes(ha(2*(iMode)));
        histfit(rmoutliers(abs(lambdaTendEnsembleDMD(iMode,:,iChem)),...
            "percentiles",[10 90]),50);
        set(gca,'LineWidth',4,'FontSize',fontSize)
        
    end
end
