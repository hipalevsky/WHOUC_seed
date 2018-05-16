function [data, B, residmean, R2] = MLR_fromGOSHIPdata_netcdf(filename)

%% Note that this is developed for AR07E bottle data provided via email from the British Oceanic Data Center
%%%% This should be similar formatting for other GO-SHIP netcdf datasets but will
%%%% need to be customized (for instance with lat/lon boundary values
%%%% prescribed below)

%%%% INPUTS:
%%%% filename = string to point to the file with all bottle data for the cruise

%%%% OUTPUTS:
%%%% data = structure with indvidual variables
%%%% B = MLR coefficients for DIC and NO3 using predictor variables: 1, T, S, O2, AOU
%%%% residmean = mean of absolute value of residuals from the MLR (for DIC and NO3)
%%%% R2 = R2 statistic for the MLR (DIC and NO3)

%% Read in data, remove flagged points, and pull out individual variables

%Pull out individual parameters
data.datetime = ncread(filename,'time')'; %seconds --> need to convert into yyyymmdd and time separately if want to be comparable
data.stn = ncread(filename,'statnum')'; %station number
data.lat = ncread(filename,'lat')'; %deg N
data.lon = ncread(filename,'lon')'; %deg E
data.pres = ncread(filename,'upress')'; %dbar
data.T = ncread(filename,'upotemp')'; %deg C
data.S = ncread(filename,'upsal')'; %pss-78
data.O2 = ncread(filename,'botoxy')'; %umol/kg
    data.botoxy_flag = ncread(filename,'botoxyflag')'; %2-9 flag
    ind = find(data.botoxy_flag > 2 & data.botoxy_flag < 5);
    data.O2(ind) = NaN;
data.NO3 = ncread(filename,'totnit')'; %umol/l
    data.totnit_flag = ncread(filename,'totnit_flag')'; %2-9 flag
    ind = find(data.totnit_flag > 2 & data.totnit_flag < 5);
    data.NO3(ind) = NaN;
data.Si = ncread(filename,'silc')'; %umol/l
    data.silc_flag = ncread(filename,'silc_flag')'; %2-9 flag
    ind = find(data.silc_flag > 2 & data.silc_flag < 5);
    data.Si(ind) = NaN;
data.PO4 = ncread(filename,'phos')'; %umol/l
    data.phos_flag = ncread(filename,'phos_flag')'; %2-9 flag
    ind = find(data.phos_flag > 2 & data.phos_flag < 5);
    data.PO4(ind) = NaN;
data.TA = ncread(filename,'alk')'; %umol/kg
    data.TA_flag = ncread(filename,'alk_flag')'; %2-9 flag
    ind = find(data.TA_flag > 2 & data.TA_flag < 5);
    data.TA(ind) = NaN;
data.DIC = ncread(filename,'dic')'; %umol/kg
    data.DIC_flag = ncread(filename,'dic_flag')'; %2-9 flag
    ind = find(data.DIC_flag > 2 & data.DIC_flag < 5);
    data.DIC(ind) = NaN;
    
%% Calculate AOU
O2equil = gsw_O2sol_SP_pt(data.S, data.T);
data.AOU = O2equil - data.O2;

%% Plot sections
%Remove station 2 (test station) at 51 N (keep data north of 55) - mostly
%E-W section so can plot by lon
    ind = find(data.lat > 55);
    M = 200;
    ymin = 0; ymax = 3500;

figure(10); clf
    subplot(321)
scatter(data.lon(ind),data.pres(ind),M,data.T(ind),'.'); colorbar;
set(gca, 'YDir', 'reverse'); ylim([ymin ymax])
title('Temperature')
    subplot(322)
scatter(data.lon(ind),data.pres(ind),M,data.S(ind),'.'); colorbar;
set(gca, 'YDir', 'reverse'); ylim([ymin ymax])
title('Salinity')
    subplot(323)
scatter(data.lon(ind),data.pres(ind),M,data.O2(ind),'.'); colorbar;
set(gca, 'YDir', 'reverse'); ylim([ymin ymax])
title('Oxygen')
    subplot(324)
scatter(data.lon(ind),data.pres(ind),M,data.AOU(ind),'.'); colorbar;
set(gca, 'YDir', 'reverse'); ylim([ymin ymax])
title('AOU')
    subplot(325)
scatter(data.lon(ind),data.pres(ind),M,data.NO3(ind),'.'); colorbar;
set(gca, 'YDir', 'reverse'); ylim([ymin ymax])
title('Nitrate')
    subplot(326)
scatter(data.lon(ind),data.pres(ind),M,data.DIC(ind),'.'); colorbar;
set(gca, 'YDir', 'reverse'); ylim([ymin ymax])
title('DIC')


%% Look at individual variable-variable relationships

figure(11); clf
    subplot(241)
scatter(data.T(ind), data.S(ind), M, data.pres(ind),'.'); colorbar; ylim([34.5 35.6])
xlabel('T'); ylabel('S');
    subplot(242)
scatter(data.T(ind), data.AOU(ind), M, data.pres(ind),'.'); colorbar
xlabel('T'); ylabel('AOU');
    subplot(243)
scatter(data.T(ind), data.NO3(ind), M, data.pres(ind),'.'); colorbar
xlabel('T'); ylabel('Nitrate');
    subplot(244)
scatter(data.T(ind), data.DIC(ind), M, data.pres(ind),'.'); colorbar
xlabel('T'); ylabel('DIC');
    subplot(245)
scatter(data.DIC(ind), data.O2(ind), M, data.pres(ind),'.'); colorbar
xlabel('DIC'); ylabel('Oxygen');
    subplot(246)
scatter(data.DIC(ind), data.AOU(ind), M, data.pres(ind),'.'); colorbar
xlabel('DIC'); ylabel('AOU');
    subplot(247)
scatter(data.NO3(ind), data.AOU(ind), M, data.pres(ind),'.'); colorbar
xlabel('Nitrate'); ylabel('AOU');
    subplot(248)
scatter(data.DIC(ind), data.NO3(ind), M, data.pres(ind),'.'); colorbar
xlabel('DIC'); ylabel('Nitrate');

%%% Note that the strongest relationships are AOU:DIC, AOU:NO3 and NO3:DIC

%% Calculate multiple linear regression (predictor variables are T, S, O2, and AOU)

%%% Remove surface-most samples (< 30 m) and those from test cast (not
%%% along section) and those from westernmost edge (freshwater)
    ind = find(data.pres > 30 & data.lat > 55 & data.lon > -42.5);

[B(:,1),BINT,R,RINT,STATS] = regress(data.DIC(ind),[ones(size(data.T(ind))) data.T(ind) data.S(ind) data.O2(ind) data.AOU(ind)]);
residmean(1) = nanmean(abs(R));
R2(1) = STATS(1);

figure(12); clf
    subplot(211)
scatter(data.pres(ind), R, [], data.lon(ind), '.'); c = colorbar;
xlabel('Pressure'); ylabel('Residuals'); title('MLR for DIC'); xlabel(c,'Longitude')

[B(:,2),BINT,R,RINT,STATS] = regress(data.NO3(ind),[ones(size(data.T(ind))) data.T(ind) data.S(ind) data.O2(ind) data.AOU(ind)]);
residmean(2) = nanmean(abs(R));
R2(2) = STATS(1);
	subplot(212)
scatter(data.pres(ind), R, [], data.lon(ind), '.'); c = colorbar;
xlabel('Pressure'); ylabel('Residuals'); title('MLR for Nitrate'); xlabel(c,'Longitude')

