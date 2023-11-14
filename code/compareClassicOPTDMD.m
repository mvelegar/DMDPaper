% DMD for 2 months of preprocessed data: 
%   (i) reconstruct the time series using different DMD
%   algorithms. Results for Classic DMD, OPTDMD with no constraints, 
%   OPTDMD with eigen values constrained to the imaginary axis, and  
%   OPTDMD with eigen values constrained to the left half plane 
%   are computed. 
%   (ii) save results for plotting
% REQUIRES:
%   (i) optdmd scripts
%   (ii) data folder with preprocessed data files

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

% Time info - for the first 3 months of data
months=201307:201309;
% Discard the first day or the first 71 snapshots
nSnapsDay = 72; % For snapshots every 20-min
% Month 1 has 30 days, 2 31 days and 3 30 days
nSnaps=load('../data/nSnaps.mat'); nSnaps=nSnaps.nSnaps;
nSnapsTotal = sum(nSnaps); nDays=nSnapsTotal/nSnapsDay;

% The chemical species info
% The 6 chemical species of interest
chem_species=cellstr(...
    ['NO  ';
    'O3  ';
    'NO2 ';
    'OH  ';
    'ISOP';
    'CO  ';]);

%% Load the preprocessed data
% Plots for deciding how to preprocess the data etc
% Figure (2) data saved here as well
% The data matrices for saving the data for preprocessing plot
% Saving 3 months of data for the 12 lats,
% with steps in this order:
% Shifted Raw data 
AStartShift3Months = load('../data/AStartShift3Months.mat');
AStartShift3Months = AStartShift3Months.AStartShift3Months;

ATendShift3Months = load('../data/ATendShift3Months.mat');
ATendShift3Months = ATendShift3Months.ATendShift3Months;

%% Comparison of DMD and OPTDMD
% fixed lev=1, fixed lat=30. 
clc; close all;
nTrainDays=40;

% Decided on these number of modes based on observations
rStart=[25 25 25 25 25 50];
rTend=[20 20 50 50 20 20]; iLat=12;

%% The optdmd setup
addpath('optdmd-master/src/');

%% The DMD algorithms

for iChem=1:length(chem_species)

    % Original data: time x nLon x nlat x chemSpecies
    YStart=AStartShift3Months(:,:,iLat,iChem);
    YTend=ATendShift3Months(:,:,iLat,iChem);


    % Now Dim 1: Lon, Dim 2: time
    YStart=permute(YStart,[2 1]);
    YTend=permute(YTend,[2 1]);

    save(['../results/YStart',chem_species{iChem},'.mat'],...
        'YStart','-v7.3');
    save(['../results/YTend',chem_species{iChem},'.mat'],...
        'YTend','-v7.3');


    % Time vector
    t=linspace(0,nTrainDays,nSnapsDay*nTrainDays); % Recon+Pred
    % sample time
    nSnapsTrain=nSnapsDay*nTrainDays; 
    % Take only the training data
    YStartTrain=YStart(:,1:nSnapsTrain);
    YTendTrain=YTend(:,1:nSnapsTrain);

    save(['../results/YStartTrain',chem_species{iChem},'.mat'],...
        'YStartTrain','-v7.3');
    save(['../results/YTendTrain',chem_species{iChem},'.mat'],...
        'YTendTrain','-v7.3');

    %% START Data DMD
    % Classic DMD
    X1=YStartTrain(:,1:nSnapsTrain-1); X2=YStartTrain(:,2:nSnapsTrain);
    % Modes
    [U,S,V]=svd(X1,'econ');
    % Reduced order
    Ur=U(:,1:rStart(iChem)); Sr=S(1:rStart(iChem),1:rStart(iChem));
    Vr=V(:,1:rStart(iChem));
    % DMD Modes
    Atilde=Ur'*X2*Vr/Sr;
    [Wr,D]=eig(Atilde);
    Phi=X2*Vr/Sr*Wr;

    % Reconstruction with classic DMD START

    % Amplitudes, project Snapshot 1
    x1=X1(:,1); b=Phi\x1;
    % Eigen values
    lambda=diag(D); omega=log(lambda);
    % Time dynamics
    time_dyn=zeros(rStart(iChem),nSnapsTrain-1);
    for j=1:nSnapsTrain
        time_dyn(:,j)=b.*exp(omega*t(j));
    end

    % Reconstruction
    Y1=Phi*time_dyn; Y1Start=real(Y1);
    save(['../results/omegaStart',chem_species{iChem},'Lat30.mat'],...
        'omega','-v7.3');


    %% OPTDMD
    % 1 --- fit to unprojected data
    imode = 1;

    opts=varpro_opts('maxiter',300);
    [w,e,b] = optdmd(YStartTrain,t(1:nSnapsTrain),rStart(iChem),imode,opts);
    clc;
    save(['../results/wStart',chem_species{iChem},'Lat30.mat'],'w','-v7.3');
    save(['../results/eStart',chem_species{iChem},'Lat30.mat'],'e','-v7.3');
    save(['../results/bStart',chem_species{iChem},'Lat30.mat'],'b','-v7.3');

    %%
    % reconstructed values
    Y2 = w*diag(b)*exp(e*t(1:nSnapsTrain)); Y2Start=real(Y2);


    %% OPTDMD - with Linear constraints
    % 1 --- fit to unprojected data
    imode = 1;
    % bounds

    % note that the first ia bounds apply to the real part of
    % alpha and the second ia bounds apply to the imaginary part
    % For unbounded, use the appropriate choice of +/- Inf

    % the below has the effect of constraining the alphas to the
    % imaginary axis

    lbc = [zeros(size(e)); -Inf*ones(size(e))];
    ubc = [zeros(size(e)); Inf*ones(size(e))];
    copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);

    opts=varpro_opts('maxiter',300);
    [w,e,b] = optdmd(YStartTrain,t(1:nSnapsTrain),rStart(iChem),imode,opts,[],[],copts);
    % reconstructed values
    Y3 = w*diag(b)*exp(e*t(1:nSnapsTrain)); Y3Start=real(Y3);
    clc;
    save(['../results/wStartCon1',chem_species{iChem},'Lat30.mat'],'w','-v7.3');
    save(['../results/eStartCon1',chem_species{iChem},'Lat30.mat'],'e','-v7.3');
    save(['../results/bStartCon1',chem_species{iChem},'Lat30.mat'],'b','-v7.3');

    % the below has the effect of constraining the alphas to the
    % left have plane
    lbc = [-Inf*ones(size(e)); -Inf*ones(size(e))];
    ubc = [zeros(size(e)); Inf*ones(size(e))];

    copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);

    opts=varpro_opts('maxiter',300);
    [w,e,b] = optdmd(YStartTrain,t(1:nSnapsTrain),rStart(iChem),imode,opts,[],[],copts);
    % reconstructed values
    Y4 = w*diag(b)*exp(e*t(1:nSnapsTrain)); Y3Start=real(Y4);
    clc;
    save(['../results/wStartCon2',chem_species{iChem},'Lat30.mat'],'w','-v7.3');
    save(['../results/eStartCon2',chem_species{iChem},'Lat30.mat'],'e','-v7.3');
    save(['../results/bStartCon2',chem_species{iChem},'Lat30.mat'],'b','-v7.3');

    %%
    save(['../results/Y1Start',chem_species{iChem},'.mat'],'Y1','-v7.3');
    save(['../results/Y2Start',chem_species{iChem},'.mat'],'Y2','-v7.3');
    save(['../results/Y3Start',chem_species{iChem},'.mat'],'Y3','-v7.3');
    save(['../results/Y4Start',chem_species{iChem},'.mat'],'Y4','-v7.3');
    clear X1 X2 Y1 Y2 Y3 Y4 Phi lambda omega


    %% TEND Data DMD

    % Classic DMD
    X1=YTendTrain(:,1:nSnapsTrain-1); X2=YTendTrain(:,2:nSnapsTrain);
    % Modes
    [U,S,V]=svd(X1,'econ');

    % Reduced order
    Ur=U(:,1:rTend(iChem)); Sr=S(1:rTend(iChem),1:rTend(iChem));
    Vr=V(:,1:rTend(iChem));
    % DMD Modes
    Atilde=Ur'*X2*Vr/Sr;
    [Wr,D]=eig(Atilde);
    Phi=X2*Vr/Sr*Wr;

    %%

    % Amplitudes, project Snapshot 1
    x1=X1(:,1); b=Phi\x1;
    % Eigen values
    lambda=diag(D); omega=log(lambda);
    % Time dynamics
    time_dyn=zeros(rTend(iChem),nSnapsTrain-1);
    for j=1:nSnapsTrain
        time_dyn(:,j)=b.*exp(omega*t(j));
    end

    % Reconstruction
    Y1=Phi*time_dyn; Y1Tend=real(Y1);
    save(['../results/omegaTend',chem_species{iChem},'Lat30.mat'],...
        'omega','-v7.3');


    %% OPTDMD
    % 1 --- fit to unprojected data
    imode = 1;
    opts=varpro_opts('maxiter',300);
    [w,e,b] = optdmd(YTendTrain,t(1:nSnapsTrain),rTend(iChem),imode,opts);
    clc;
    save(['../results/wTend',chem_species{iChem},'Lat30.mat'],'w','-v7.3');
    save(['../results/eTend',chem_species{iChem},'Lat30.mat'],'e','-v7.3');
    save(['../results/bTend',chem_species{iChem},'Lat30.mat'],'b','-v7.3');
    % reconstructed values
    Y2 = w*diag(b)*exp(e*t(1:nSnapsTrain)); Y2Tend=real(Y2);
    %% OPTDMD with Linear Constraints

    lbc = [zeros(size(e)); -Inf*ones(size(e))];
    ubc = [zeros(size(e)); Inf*ones(size(e))];

    copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);


    opts=varpro_opts('maxiter',300);
    [w,e,b] = optdmd(YTendTrain,t(1:nSnapsTrain),rTend(iChem),imode,opts,[],[],copts);
    clc;
    % reconstructed values
    Y3 = w*diag(b)*exp(e*t(1:nSnapsTrain)); Y3Tend=real(Y3);
    save(['../results/wTendCon1',chem_species{iChem},'Lat30.mat'],'w','-v7.3');
    save(['../results/eTendCon1',chem_species{iChem},'Lat30.mat'],'e','-v7.3');
    save(['../results/bTendCon1',chem_species{iChem},'Lat30.mat'],'b','-v7.3');

    lbc = [-Inf*ones(size(e)); -Inf*ones(size(e))];
    ubc = [zeros(size(e)); Inf*ones(size(e))];

    copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);

    opts=varpro_opts('maxiter',300);
    [w,e,b] = optdmd(YTendTrain,t(1:nSnapsTrain),rTend(iChem),imode,opts,[],[],copts);
    clc;
    % reconstructed values
    Y4 = w*diag(b)*exp(e*t(1:nSnapsTrain)); Y4Tend=real(Y4);
    save(['../results/wTendCon2',chem_species{iChem},'Lat30.mat'],'w','-v7.3');
    save(['../results/eTendCon2',chem_species{iChem},'Lat30.mat'],'e','-v7.3');
    save(['../results/bTendCon2',chem_species{iChem},'Lat30.mat'],'b','-v7.3');

    save(['../results/Y1Tend',chem_species{iChem},'.mat'],'Y1','-v7.3');
    save(['../results/Y2Tend',chem_species{iChem},'.mat'],'Y2','-v7.3');
    save(['../results/Y3Tend',chem_species{iChem},'.mat'],'Y3','-v7.3');
    save(['../results/Y4Tend',chem_species{iChem},'.mat'],'Y4','-v7.3');




end
clc;
