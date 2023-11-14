% Plot the spatial modes for all latitudes 
% Requires: computeSpatialModes results from results folder


%% Plot spatial modes for all latitudes for both Tend and Start data
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
% Isolate the land cells from ocean cells
landMapLim = landmap(latVecIndLim(1):latVecIndLim(2),:);

% Time info
nSnapsDay = 72; % For snapshots every 20-min
nTrainDays=40; nSnapsTrain=nTrainDays*nSnapsDay;
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
%% No of Modes 
rStart=[25 25 25 25 25 50]; 
rTend=[20 20 50 50 20 20]; 
% For the plots
addpath('cmocean_v1.4/cmocean/');

%% Load the spatial modes
% For OPT DMD with Linear Constraints (left plane evals only)
wStartCon2OPTAll=load('../results/wStartOPTAll.mat');
wStartCon2OPTAll=wStartCon2OPTAll.wStartCon2OPTAll; 
wStartCon2OPTAll=real(wStartCon2OPTAll);
wTendCon2OPTAll=load('../results/wTendOPTAll.mat');
wTendCon2OPTAll=wTendCon2OPTAll.wTendCon2OPTAll;
wTendCon2OPTAll=real(wTendCon2OPTAll);

nModes=10; 
wStartOPTPlt=NaN(nLon,nlat,nModes,length(chem_species)); 
wTendOPTPlt=NaN(nLon,nlat,nModes,length(chem_species));
wStartCon2OPTPlt=NaN(nLon,nlat,nModes,length(chem_species)); 
wTendCon2OPTPlt=NaN(nLon,nlat,nModes,length(chem_species));

for iChem=1:length(chem_species)
   for iMode=1:nModes
       temp=wStartCon2OPTAll(:,iMode,iChem,:); temp=squeeze(temp);
       wStartCon2OPTPlt(:,:,iMode,iChem)=temp;
       temp=wTendCon2OPTAll(:,iMode,iChem,:); temp=squeeze(temp);
       wTendCon2OPTPlt(:,:,iMode,iChem)=temp;
   end
end

%% Plot the spatial modes
 cmap = cmocean('balance'); fontSize=18;
 [X2d,Y2d]=meshgrid(x,yLim);
for iChem=1:length(chem_species)
    figure();
    ha = tight_subplot(4,2,[.01 .09],[.1 0.01],[.05 .09]);
    for iMode=1:4
       % START PLOTS
       cmin=min(min(min(wStartCon2OPTPlt(:,:,1:4,iChem)))); 
       cmax=max(max(max(wStartCon2OPTPlt(:,:,1:4,iChem))));
       axes(ha(2*(iMode-1)+1));
       pcolor(X2d,Y2d,wStartCon2OPTPlt(:,:,2*(iMode-1)+1,iChem)'); grid off; 
       shading interp; hold on;
       [~,h]=contour(X2d,Y2d,landMapLim,[0.9999 0.9999]);h.LineColor='k'; 
       h.LineWidth=2; hold off;
       colormap(cmap); 
       % center the colormap
       cLimCenter=max(abs([cmin cmax]));
       caxis([-cLimCenter cLimCenter]);
       hp1 = get(ha(1),'Position'); hp7 = get(ha(7),'Position');
       cb_x = hp7(1) + hp7(3) + 0.005;
       cb_y = hp7(2);
       cb_w = 0.025;
       cb_h = hp1(2) + hp1(4) - hp7(2);
       colorbar('Position', [cb_x cb_y cb_w cb_h]);
       
        if (2*(iMode-1)+1)==1 || (2*(iMode-1)+1)==3 || (2*(iMode-1)+1)==5 ...
                || (2*(iMode-1)+1)==7
            set(gca,'LineWidth',4,...
                'ytick',[-14,6,26],'YTickLabel',...
                ({'-14','6','26'}),'LineWidth',4,...
                'FontSize',fontSize,'XTickLabel',...
                []);
        end
        if (2*(iMode-1)+1)==7 
            set(gca,...
                'xtick',[-90,0,90],'LineWidth',4,...
                'FontSize',fontSize,'XTickLabel',...
                ({'-90','0','90'}));
       end
        if (2*(iMode-1)+1)==1 || (2*(iMode-1)+1)==3 || (2*(iMode-1)+1)==5
            xlabel('$\mathrm{Lon}$','Interpreter','Latex','Fontsize',fontSize);
        end
       
       % TEND PLOTS
       cmin=min(min(min(wTendCon2OPTPlt(:,:,1:4,iChem)))); 
       cmax=max(max(max(wTendCon2OPTPlt(:,:,1:4,iChem))));
       axes(ha(2*iMode));
       pcolor(X2d,Y2d,wTendCon2OPTPlt(:,:,2*(iMode-1)+1,iChem)'); grid off; 
       shading interp; hold on;
       [~,h]=contour(X2d,Y2d,landMapLim,[0.9999 0.9999]);h.LineColor='k'; 
       h.LineWidth=2; hold off;
       colormap(cmap); 
       % center the colormap
       cLimCenter=max(abs([cmin cmax]));
       caxis([-cLimCenter cLimCenter]);
       hp1 = get(ha(1),'Position'); hp8 = get(ha(8),'Position');
       cb_x = hp8(1) + hp8(3) + 0.005;
       cb_y = hp8(2);
       cb_w = 0.025;
       cb_h = hp1(2) + hp1(4) - hp8(2);
       colorbar('Position', [cb_x cb_y cb_w cb_h]);
       if  (2*iMode)==8
            set(gca,...
                'xtick',[-90,0,90],'LineWidth',4,...
                'FontSize',fontSize,'XTickLabel',...
                ({'-90','0','90'}));
       end
       
       if (2*iMode)==2||(2*iMode)==4||(2*iMode)==6||(2*iMode)==8
          set(gca,...
                'LineWidth',4,...
                'FontSize',fontSize,'YTickLabel',...
                []); 
       end
       
       
    end
end

