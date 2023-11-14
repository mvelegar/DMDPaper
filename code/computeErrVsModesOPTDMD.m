% DMD for 40 days of preprocessed data:
%   (i) reconstruct the time series using OPT DMD
%   algorithm with eigen values constrained to the left half plane using
%   different number of modes
%   (ii) compute the relative errors, and save results for plotting
% REQUIRES:
%   (i) optdmd scripts
%   (ii) data folder with preprocessed data files

clear variables; close all; clc

%% Set up grid, time and species variables
% Landmap information
ncfileMap='../data/landmap.nc';
% Read the Spatial grid
x=ncread(ncfileMap,'lon');% Longitude(-180:5:175)

% Time info - for the first 40 days of data
% Time vector
nSnapsDay=72;
nTrainDays=40; nSnapsTrain=nSnapsDay*nTrainDays;
t=linspace(0,nTrainDays,nSnapsTrain);


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
YStartTrainAll = NaN(length(x),nSnapsTrain,length(chem_species));
YTendTrainAll = NaN(length(x),nSnapsTrain,length(chem_species));
for iChem = 1:length(chem_species)
    % Training Data
    Y = load(['../results/YStartTrain',chem_species{iChem},'.mat']);
    Y = Y.YStartTrain; YStartTrainAll(:,:,iChem)=real(Y); clear Y;
    Y = load(['../results/YTendTrain',chem_species{iChem},'.mat']);
    Y = Y.YTendTrain; YTendTrainAll(:,:,iChem)=real(Y);
end

%% The optdmd setup
addpath('optdmd-master/src/');


%% The relative errors
clc; close all;
R=1:50;

relErrStart=NaN(length(R),length(chem_species));
relErrTend=NaN(length(R),length(chem_species));
for iChem=1:length(chem_species)
    relErrStart=NaN(length(R),1);
    YStart=YStartTrainAll(:,:,iChem);
    YTend=YTendTrainAll(:,:,iChem);
    for iRank=1:length(R)

        %         OPT-DMD
        %         target rank
        r = R(iRank);
        % OPTDMD with Linear Constraints
        lbc = [-Inf*ones(r,1); -Inf*ones(r,1)];
        ubc = [zeros(r,1); Inf*ones(r,1)];
        imode = 1;
        copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);
        opts=varpro_opts('maxiter',100);


        % Relative Error for START
        [w,e,b] = optdmd(YStart,t(1:nSnapsTrain),r,imode,opts,[],[],copts);
        clc;
        % Reconstruct
        Y2 = w*diag(b)*exp(e*t); Y2=real(Y2);
        %        Relative 2 norm error
        relErrStart(iRank,iChem)=norm(YStart-Y2)/norm(YStart);
        clear w e b Y2;

        % Relative error for TEND
        [w,e,b] = optdmd(YTend,t(1:nSnapsTrain),r,imode,opts,[],[],copts);
        clc;
        % Reconstruct
        Y2 = w*diag(b)*exp(e*t); Y2=real(Y2);
        %        Relative 2 norm error
        relErrTend(iRank,iChem)=norm(YTend-Y2)/norm(YTend);
        clear w e b Y2;

        clc;
    end
end
save(['../results/relErrVsRStartAllLat30.mat'],'relErrStart','-v7.3');
save(['../results/relErrVsRTendAllLat30.mat'],'relErrTend','-v7.3');

