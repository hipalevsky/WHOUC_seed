function [data, B, residmean, R2] = MLR_fromGOSHIPdata(filename, x)

%% Note that this is developed for AR07E bottle data downloaded from CLIVAR database (https://www.nodc.noaa.gov/ocads/oceans/RepeatSections/clivar_data.html)
%%%% This should be similar formatting for other CLIVAR datasets but will
%%%% need to be customized (for instance with lat/lon boundary values
%%%% prescribed below)

%%%% INPUTS:
%%%% filename = string to point to the file with all bottle data for the cruise
%%%% x = factor to adjust for different column arrangements

%%%% OUTPUTS:
%%%% data = structure with indvidual variables
%%%% B = MLR coefficients for DIC and NO3 using predictor variables: 1, T, S, O2, AOU
%%%% residmean = mean of absolute value of residuals from the MLR (for DIC and NO3)
%%%% R2 = R2 statistic for the MLR (DIC and NO3)

%% Read in data, remove filler/flagged points, and pull out individual variables
datain = xlsread(filename);

%Remove points with -999 as filler data and replace with NaN
[m,n] = size(datain);
    for i = 1:n
        datain((datain(:,i) == -999),i) = NaN;
        if i > 3 %Remove points with questionable data (3 flag)
            datain((datain(:,i) == 3),i-1) = NaN;
        end
            A = find(datain(:,i) == 3);
            num3s(i) = length(A); %number of points removed
    end

%Pull out individual parameters
data.stn = datain(:,1);
data.date = datain(:,5); %seems like yyyymmdd, but year is 2005 and should be 2007...other things in data file also out of date, but think data are ok...
data.time = datain(:,6); 
data.lat = datain(:,7); 
data.lon = datain(:,8); 
data.pres = datain(:,10); %dbar, from CTD
data.T = datain(:,11); %ITS-90, from CTD
data.S = datain(:,12); %PSS-78, from CTD
data.O2 = datain(:,16 - x); %umol/kg, from bottle
data.Si = datain(:,18 - x); %umol/kg, from bottle
data.NO3 = datain(:,20 - x); %umol/kg, from bottle
data.PO4 = datain(:,22); %umol/kg, from bottle
data.DIC = datain(:,24); %umol/kg, from bottle
data.TA = datain(:,26); %umol/kg, from bottle

%% Calculate AOU
O2equil = gsw_O2sol_SP_pt(data.S, data.T);
data.AOU = O2equil - data.O2;

%% Plot sections
%Remove station 2 (test station) at 51 N (keep data north of 55) - mostly
%E-W section so can plot by lon
    ind = find(data.lat > 55);
    M = 200;
    ymin = 0; ymax = 3500;

figure(1); clf
    subplot(321)
scatter(data.lon(ind),data.pres(ind),M,data.T(ind),'.'); colorbar;
set(gca, 'YDir', 'reverse'); ylim([ymin ymax])
title('Temperature (ITS-90)')
xlabel('Longitude'); ylabel('Depth (m)')
    subplot(322)
scatter(data.lon(ind),data.pres(ind),M,data.S(ind),'.'); colorbar;
set(gca, 'YDir', 'reverse'); ylim([ymin ymax])
title('Salinity (PSS-78)'); caxis([34.5 35.5])
xlabel('Longitude'); ylabel('Depth (m)')
    subplot(323)
scatter(data.lon(ind),data.pres(ind),M,data.O2(ind),'.'); colorbar;
set(gca, 'YDir', 'reverse'); ylim([ymin ymax])
title('Oxygen (\mumol/kg)')
xlabel('Longitude'); ylabel('Depth (m)')
    subplot(324)
scatter(data.lon(ind),data.pres(ind),M,data.AOU(ind),'.'); colorbar;
set(gca, 'YDir', 'reverse'); ylim([ymin ymax])
title('AOU (\mumol/kg)')
xlabel('Longitude'); ylabel('Depth (m)')
    subplot(325)
scatter(data.lon(ind),data.pres(ind),M,data.NO3(ind),'.'); colorbar;
set(gca, 'YDir', 'reverse'); ylim([ymin ymax])
title('Nitrate (\mumol/kg)')
xlabel('Longitude'); ylabel('Depth (m)')
    subplot(326)
scatter(data.lon(ind),data.pres(ind),M,data.DIC(ind),'.'); colorbar;
set(gca, 'YDir', 'reverse'); ylim([ymin ymax])
title('DIC (\mumol/kg)')
xlabel('Longitude'); ylabel('Depth (m)')


%% Look at individual variable-variable relationships

figure(2); clf
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

figure(3); clf
    subplot(211)
scatter(data.pres(ind), R, [], data.lon(ind), '.'); c = colorbar;
xlabel('Pressure'); ylabel('Residuals'); title('MLR for DIC'); xlabel(c,'Longitude')

[B(:,2),BINT,R,RINT,STATS] = regress(data.NO3(ind),[ones(size(data.T(ind))) data.T(ind) data.S(ind) data.O2(ind) data.AOU(ind)]);
residmean(2) = nanmean(abs(R));
R2(2) = STATS(1);
	subplot(212)
scatter(data.pres(ind), R, [], data.lon(ind), '.'); c = colorbar;
xlabel('Pressure'); ylabel('Residuals'); title('MLR for Nitrate'); xlabel(c,'Longitude')

