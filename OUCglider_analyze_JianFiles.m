%% Load gridded data shared by Jian
%Note that if he has comparable data structures for other missions, this
%should be easy to scale up
load C:/Users/Hilary/Dropbox/WHOUCSeedFundingProject/GliderData/Mission6_WHOUC3/WHOUC3_profiledata_201709_201712_jian.mat

%% Plot lat-lon map over time
profiledata.time = datenum(profiledata.yyddmm);
profiledata.latlon_plot = profiledata.latlon/100;

figure(1); clf
scatter(profiledata.latlon_plot(:,2), profiledata.latlon_plot(:,1), [], (profiledata.time - min(profiledata.time)), 'filled')
h = colorbar; ylabel(h, 'Days into mission')
xlabel('Longitude'); ylabel('Latitude'); title('Map of glider profile locations, Mission 6')

%% Calculate derivative variables (gibbs seawater toolbox)
[profiledata.standSA, ~] = gsw_SA_from_SP(profiledata.standsalt, profiledata.standP, profiledata.latlon_plot(:,2), profiledata.latlon_plot(:,1)); %absolute salinity
profiledata.standConservT = gsw_CT_from_t(profiledata.standSA, profiledata.standT, profiledata.standP); %conservative temperature
profiledata.standSigma0 = gsw_sigma0(profiledata.standSA,profiledata.standConservT); %potential density anomaly at surface
profiledata.standO2equil = gsw_O2sol(profiledata.standSA, profiledata.standConservT, profiledata.standP, profiledata.latlon_plot(:,2), profiledata.latlon_plot(:,1)); %O2 at equilibrium, surf ocean

%% Calculate raw O2 derivative variables and plot time series sections
%For now this is just a unit conversion of raw data - needs gain and drift correction (below)
profiledata.standO2_raw = profiledata.standO2./((profiledata.standSigma0 + 1000)/1000); %umol/L to umol/kg, normalized to density at surface (sigma0)
profiledata.standAOU_raw = profiledata.standO2equil - profiledata.standO2_raw; %apparent oxygen utilization (umol/kg)
profiledata.standO2biosat_raw = (profiledata.standO2_raw./profiledata.standO2equil - 1)*100; %O2 % biological supersaturation

figure(2); clf
    subplot(511)
imagesc(profiledata.standConservT); colorbar;
title('Temperature, deg C');
    subplot(512)
imagesc(profiledata.standSA); colorbar;
title('Salinity');
    subplot(513)
imagesc(profiledata.standO2_raw); colorbar;
title('Uncorr. O_2 concentration (\mumol/kg)');
    subplot(514)
imagesc(profiledata.standAOU_raw); colorbar;
title('Uncorr. AOU (\mumol/kg)');
    subplot(515)
imagesc(profiledata.standO2biosat_raw); colorbar;
title('Uncorr. O_2 biosat (%)');

%% Illustrate the offset in the oxygen data by plotting surface-most points
%Quick and dirty assumption would be that AOU and O2_biosat both = 0 at
%surface (i.e. in equilibrium with atmosphere) --> This would preclude
%interpretation of surface data, but is a rough first step
%This is supported roughly by AR07E occupation from Sept. 2007, for which
%mean AOU of top 20 m is -2.8 +/- 6.6 umol/kg

%Calculate linear regression to correct AOU to be zero at surface
%throughout
P = polyfit(profiledata.time - min(profiledata.time), profiledata.standAOU_raw(1,:)', 1);

figure(3); clf
    subplot(211)
plot(profiledata.time, profiledata.standAOU_raw(1,:), 'k.'); hold on;
plot(profiledata.time, P(2) + (profiledata.time - min(profiledata.time)).*P(1), 'r--'); hold on;
datetick('x'); title('Surface uncorr. AOU over mission')
    subplot(212)
plot(profiledata.time, profiledata.standO2biosat_raw(1,:), 'k.');
datetick('x'); title('Surface uncorr. O_2 biosat over mission')

%% Determine gain and offset corrections, and use to correct O2 data
gain = (profiledata.standO2_raw(1,1) + P(2))/profiledata.standO2_raw(1,1); %based on surface O2 conc from 1st profile and intercept of AOU linear fit
drift = P(1); %umol/kg/day

    [numdepths, ~] = size(profiledata.standO2_raw);

profiledata.standO2_corr = profiledata.standO2_raw.*gain + drift.*repmat((profiledata.time - min(profiledata.time)), 1, numdepths)'; %correct for gain and drift
profiledata.standAOU_corr = profiledata.standO2equil - profiledata.standO2_corr; %apparent oxygen utilization (umol/kg)
profiledata.standO2biosat_corr = (profiledata.standO2_corr./profiledata.standO2equil - 1)*100; %O2 % biological supersaturation

figure(4); clf
    subplot(511)
imagesc(profiledata.standConservT); colorbar;
title('Conservative temperature');
ylabel('dbar')
    subplot(512)
imagesc(profiledata.standSA); colorbar;
title('Absolute salinity');
ylabel('dbar')
    subplot(513)
imagesc(profiledata.standO2_corr); colorbar;
title('Corrected O_2 concentration (\mumol/kg)');
ylabel('dbar')
    subplot(514)
imagesc(profiledata.standAOU_corr); colorbar;
ylabel('dbar')
title('Corrected AOU (\mumol/kg)');
    subplot(515)
imagesc(profiledata.standO2biosat_corr); colorbar;
title('Corrected O_2 biosat (%)');
ylabel('dbar')
xlabel('Profile number')

%% Show how O2 was corrected to surface value at equilibrium

figure(5); clf
    subplot(211)
histogram(profiledata.standAOU_corr(1,:)); hold on;
title('Histogram of corrected surface AOU over mission')
xlim([-10 10])
    subplot(212)
histogram(profiledata.standO2biosat_corr(1,:)); hold on;
title('Histogram of corrected surface O_2 biosat over mission')
xlim([-4 4])

%% Read in AR07E bottle data and calculate MLRs for nitrate and DIC
AR07E_bottledata

%% Plot location of glider mission over map of AR07E occupations
figure(6); clf
    latminplot = 48; latmaxplot = 68;
    lonminplot = -55; lonmaxplot = 0;
    M = 15;
m_proj('lambert','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
m_plot(data_05.lon, data_05.lat,'b.','markersize',M); hold on; %AR07E 2005
m_plot(data_07.lon, data_07.lat,'c.','markersize',M-1); hold on; %AR07E 2007
m_plot(data_14.lon(ind), data_14.lat(ind),'r.','markersize',M-2); hold on; %AR07E 2014
m_plot(profiledata.latlon_plot(:,2), profiledata.latlon_plot(:,1),'k.','markersize',M/1.5); hold on; %WHOUC Mission 6
%m_plot(-40-34.5/60, 59+56.7/60, 'm.','markersize',M*1.5); hold on; %OOI Irminger site
m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
legend('','AR07E 2005','AR07E 2007','AR07E 2014','WHOUC Mission 6 glider','location','southwest')

%% Calculate DIC and NO3 from glider data using 2007 AR07E MLR

profiledata.DIC_wMLR07 = B_07(1,1) + B_07(2,1)*profiledata.standT + B_07(3,1)*profiledata.standsalt + B_07(4,1)*profiledata.standO2_corr + B_07(5,1)*profiledata.standAOU_corr;
profiledata.NO3_wMLR07 = B_07(1,2) + B_07(2,2)*profiledata.standT + B_07(3,2)*profiledata.standsalt + B_07(4,2)*profiledata.standO2_corr + B_07(5,2)*profiledata.standAOU_corr;

figure(7); clf
    subplot(211)
imagesc(profiledata.DIC_wMLR07); colorbar;
ylabel('dbar')
title('DIC from AR07E-07 MLR (\mumol/kg)');
    subplot(212)
imagesc(profiledata.NO3_wMLR07); colorbar;
title('Nitrate from AR07E-07 MLR (\mumol/kg)');
ylabel('dbar')
xlabel('Profile number')


