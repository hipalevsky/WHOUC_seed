%% Extract WOA13 monthly O2 data to calculate climatological AOU

%% Read on World Ocean Atlas 2013 data to get a rough data set to work with

% Oxygen data ("all" is data from entire time period, objectively analyzed
% onto 1 degree and depth grid for every month)
for i = 1:9
ncid = netcdf.open(['woa13_all_o0' num2str(i) '_01.nc'],'NC_NOWRITE');
    lat = double(netcdf.getVar(ncid,1));
    lon = double(netcdf.getVar(ncid,3));
    depth = double(netcdf.getVar(ncid,5));
    time(i) = double(netcdf.getVar(ncid,7));
    O(:,:,:,i) = double(netcdf.getVar(ncid,9));
netcdf.close(ncid)
end

for i = 10:12
ncid = netcdf.open(['woa13_all_o' num2str(i) '_01.nc'],'NC_NOWRITE');
    time(i) = double(netcdf.getVar(ncid,7));
    O(:,:,:,i) = double(netcdf.getVar(ncid,9));
netcdf.close(ncid)
end

%Reshape to have longitudes run from 0 to 360
Odouble = [O; O];
O = Odouble(181:540,:,:,:);

% Temperature data ("A5B2" is data from 2005-2012, objectively analyzed
% onto 1 degree and depth grid for every month. Other option is "decav",
% which takes mean of six decadal means)
for i = 1:9
ncid = netcdf.open(['woa13_A5B2_t0' num2str(i) '_01.nc'],'NC_NOWRITE');
    lat = double(netcdf.getVar(ncid,1));
    lon = double(netcdf.getVar(ncid,3));
    depth = double(netcdf.getVar(ncid,5));
    time(i) = double(netcdf.getVar(ncid,7));
    T(:,:,:,i) = double(netcdf.getVar(ncid,9));
netcdf.close(ncid)
end

for i = 10:12
ncid = netcdf.open(['woa13_A5B2_t' num2str(i) '_01.nc'],'NC_NOWRITE');
    time(i) = double(netcdf.getVar(ncid,7));
    T(:,:,:,i) = double(netcdf.getVar(ncid,9));
netcdf.close(ncid)
end

%Reshape to have longitudes run from 0 to 360
Tdouble = [T; T];
T = Tdouble(181:540,:,:,:);

% Salinity data ("A5B2" is data from 2005-2012, objectively analyzed
% onto 1 degree and depth grid for every month. Other option is "decav",
% which takes mean of six decadal means)
for i = 1:9
ncid = netcdf.open(['woa13_A5B2_s0' num2str(i) '_01.nc'],'NC_NOWRITE');
    lat = double(netcdf.getVar(ncid,1));
    lon = double(netcdf.getVar(ncid,3));
    depth = double(netcdf.getVar(ncid,5));
    time(i) = double(netcdf.getVar(ncid,7));
    S(:,:,:,i) = double(netcdf.getVar(ncid,9));
netcdf.close(ncid)
end

for i = 10:12
ncid = netcdf.open(['woa13_A5B2_s' num2str(i) '_01.nc'],'NC_NOWRITE');
    time(i) = double(netcdf.getVar(ncid,7));
    S(:,:,:,i) = double(netcdf.getVar(ncid,9));
netcdf.close(ncid)
end

%Reshape to have longitudes run from 0 to 360
Sdouble = [S; S];
S = Sdouble(181:540,:,:,:);

% Revise definition of longitude to be on 0-360 rather than -180 to 180
% scale
lon = [0.5:359.5]';

%Revide definition of time to be in days
time = [15;46;75;106;136;167;197;228;259;289;320;350];

%% Calculate density and convert O2 units
[~,latgrid,depthgrid,~] = ndgrid(lon,lat,depth,time);
pres = sw_pres(depthgrid,latgrid); %pressure at the given depth and latitude, db
dens = sw_dens(S,T,pres); %density at given T, S and pressure, kg/m3
pden = sw_pden(S,T,pres,zeros(size(dens))); %potential density at given T, S and pressure, relative to the surface, kg/m3
O2 = O./[(0.0223916)*(pden./1000)]; O2 = double(O2); %convert to umol/kg from ml/L: O2 = 22391.6 ml/mol (Garcia and Gordon, 1992)

%% Pull out climatology of T, S, and O2 in Iceland Basin
    stnlat = 58.5; stnlon = -22.5+360;
woalatid = find(lat == stnlat); woalonid = find(lon == stnlon);
depthid = find(depth == 5);

O2_clim = squeeze(O2(woalonid, woalatid, depthid,:));
T_clim = squeeze(T(woalonid, woalatid, depthid,:));
S_clim = squeeze(S(woalonid, woalatid, depthid,:));

%% Calculate solubility and saturation
O2sol_clim = O2sol(S_clim,T_clim);
O2sat_clim = (O2_clim - O2sol_clim)./O2sol_clim*100;
AOU_clim = O2sol_clim - O2_clim;

%% Rough plot
M = 15; 
M2 = 5;
figure(1); clf
    subplot(211)
plot([1:12],O2sat_clim,'k.','markersize',M+M2); hold on;
plot([0:13],zeros(14,1),'k--'); hold on;
axis([0.5 12.5 -6 12]);
legend('WOA13 O_2 climatology')
ylabel('O_2 saturation (%)'); xlabel('Month')
    subplot(212)
plot([1:12],AOU_clim,'k.','markersize',M+M2); hold on;
plot([0:13],zeros(14,1),'k--'); hold on;
axis([0.5 12.5 -30 10]);
legend('WOA13 AOU climatology')
ylabel('AOU (\mumol/kg)'); xlabel('Month')

%% Create structure with data to output for location

WOAclim.O2 = squeeze(O2(woalonid, woalatid, :,:));
WOAclim.T = squeeze(T(woalonid, woalatid, :,:));
WOAclim.S = squeeze(S(woalonid, woalatid, :,:));
WOAclim.O2sat = (WOAclim.O2 - O2sol(WOAclim.S,WOAclim.T))./O2sol(WOAclim.S,WOAclim.T)*100;
WOAclim.AOU = O2sol(WOAclim.S,WOAclim.T) - WOAclim.O2;
WOAclim.depth = depth;
WOAclim.dayinyr = time;
WOAclim.latlon = [stnlat stnlon];

save WOA13clim_IcelandBasin.mat WOAclim
