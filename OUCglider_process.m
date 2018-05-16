%%% Script to extract OUC glider data, modeled off of ws_visualize for
%%% OOI Irminger Sea (pipeline for 2018 Oceanography paper)

DorR = 'c'; %added a case for OUC gliders
rootDir = 'C://Users/Hilary/Dropbox/OUC/WH_OUC3_fall2017';
[L2,meta_L2,L1,meta_L1,L0,meta_L0] = ws_process(rootDir,DorR);


%% Scatter plot of L2 data
figure(1); clf
subplot(211)
scatter(L2.daten, L2.depth_interp, [], L2.temperature);
    colormap('parula'); colorbar
    set(gca,'YDir','reverse'); datetick('x',2,'keeplimits');
    ylabel('Depth'); title('Temperature (deg C)')
subplot(212)
scatter(L2.daten, L2.depth_interp, [], L2.oxygen_concentration);
    colormap('parula'); colorbar
    set(gca,'YDir','reverse'); datetick('x',2,'keeplimits');
    ylabel('Depth'); title('Oxygen (uM, raw)')
    caxis([200 250])
    
%% Grid and visualize science data from each glider

in = L2;

%Pick which science variables to include
    scivars = [in.temperature in.salinity in.oxygen_concentration in.oxygen_saturation in.calphase_oxygen];
    var_titles = {'Temperature','Salinity','Oxygen concentration','Oxygen saturation','Oxygen calphase'};
    depth_grid = [5:10:995];

%Filter out bad data points for each variable
for j = 1:length(var_titles)
    indreal = find(real(scivars(:,j)) == scivars(:,j));   
    indbad = find(scivars(:,j) < (nanmean(scivars(indreal,j)) - 4*nanstd(scivars(indreal,j)))...
            | scivars(:,j) > (nanmean(scivars(indreal,j)) + 4*nanstd(scivars(indreal,j))));
    scivars(indbad,:) = NaN; %filter out all science data where any bad data exists (conservative choice, but doesn't seem to lose much and more effectively removes all outliers)
end

%Grid data by depth
[out] = glider_grid(in.daten,in.latitude,in.longitude,in.depth_interp,in.profile_index,in.profile_direction,scivars,depth_grid);
    out.depth_grid = depth_grid;
    out.var_titles = var_titles;
[X,Y] = meshgrid(out.time_start,depth_grid);


%Calculate MLD for plotting
pdens = sw_dens0(squeeze(out.scivars(:,2,:)),squeeze(out.scivars(:,1,:)));
    criterion = 0.03; %density change from surface
[mld03] = mld_calc(pdens,out.depth_grid,criterion);
    criterion = 0.125; %density change from surface
[mld125] = mld_calc(pdens,out.depth_grid,criterion);


%% Plot data
figure(1); clf %by depth
set(gcf,'color','w')
    x0=1;
    y0=1;
    width=19;
    height=19;
    set(gcf,'units','centimeters','position',[x0,y0,width,height])
%Find profiles with good temperature data from surface
ind = find(~isnan(squeeze(out.scivars(1,1,:))));

for i = 1:length(var_titles)
    subplot(length(var_titles),1,i)
    contourf(X(:,ind),Y(:,ind),squeeze(out.scivars(:,i,ind)),'linecolor','none'); hold on;
    plot(out.time_start(ind),mld03(ind),'color',nicecolor('kkw'),'linewidth',2); hold on;
    plot(out.time_start(ind),mld125(ind),'color',nicecolor('wwwk'),'linewidth',2); hold on;
    colormap('parula');
    colorbar
    set(gca,'YDir','reverse');
    datetick('x',2,'keeplimits');
    ylabel('Depth (m)');
    title(var_titles(i));
end
