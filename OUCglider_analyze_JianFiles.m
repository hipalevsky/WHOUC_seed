%% Read in AR07E bottle data and calculate MLRs for nitrate and DIC
AR07E_bottledata

%% Load climatological AOU from World Ocean Atlas
%Extracted using WOA_climatology_WHOUCsite.m
load WOA13clim_IcelandBasin.mat

%% Load gridded data shared by Jian

diagnostic_plots = 1; %if 1 will make diagnostic plots

for mission = [1, 2, 4, 5, 6]

if mission == 1
    load C:/Users/Hilary/Dropbox/WHOUCSeedFundingProject/GliderData/Mission1_WHOUC1/WHOUC1_profiledata_oxygen_information_201506_201511_JZ.mat
elseif mission == 2
    load C:/Users/Hilary/Dropbox/WHOUCSeedFundingProject/GliderData/Mission2_WHOUC2/WHOUC2_profiledata_oxygen_information_201512_201606_JZ.mat
elseif mission == 4
    load C:/Users/Hilary/Dropbox/WHOUCSeedFundingProject/GliderData/Mission4_WHOUC3/WHOUC3_profiledata_oxygen_information_201611_201705_JZ.mat
elseif mission == 5
    load C:/Users/Hilary/Dropbox/WHOUCSeedFundingProject/GliderData/Mission5_WHOUC2/WHOUC2_profiledata_oxygen_information_20705_201707_JZ.mat
elseif mission == 6
    load C:/Users/Hilary/Dropbox/WHOUCSeedFundingProject/GliderData/Mission6_WHOUC3/WHOUC3_profiledata_oxygen_information_201709_201712_JZ.mat
end

%% Plot lat-lon map over time
profiledata.time = datenum(profiledata.yyddmm);
%Set bad (zero) time data to NaN
    timebad = find(profiledata.time == 0);
    profiledata.time(timebad) = NaN;
profiledata.year = str2num(datestr(profiledata.time,10)); %year of sampling
profiledata.month = str2num(datestr(profiledata.time,5)); %month of sampling
profiledata.day = profiledata.time - datenum(profiledata.year,0,0); %Julian day w/o year to create annual composite
profiledata.latlon_plot = profiledata.latlon/100;

if diagnostic_plots == 1
    figure(1); clf
    scatter(profiledata.latlon_plot(:,2), profiledata.latlon_plot(:,1), [], (profiledata.time - min(profiledata.time)), 'filled')
    h = colorbar; ylabel(h, 'Days into mission')
    xlabel('Longitude'); ylabel('Latitude'); title({['Map of glider profile locations, Mission ' num2str(mission)]})
end

%% Calculate derivative variables (gibbs seawater toolbox)
[profiledata.standSA, ~] = gsw_SA_from_SP(profiledata.standsalt, profiledata.standP, profiledata.latlon_plot(:,2), profiledata.latlon_plot(:,1)); %absolute salinity
profiledata.standConservT = gsw_CT_from_t(profiledata.standSA, profiledata.standT, profiledata.standP); %conservative temperature
profiledata.standSigma0 = gsw_sigma0(profiledata.standSA,profiledata.standConservT); %potential density anomaly at surface
profiledata.standO2equil = gsw_O2sol(profiledata.standSA, profiledata.standConservT, profiledata.standP, profiledata.latlon_plot(:,2), profiledata.latlon_plot(:,1)); %O2 at equilibrium, surf ocean

%% Calculate raw O2 derivative variables and plot time series sections
%For now this is just a unit conversion of raw data - needs gain and drift correction (below)
profiledata.standO2_raw = profiledata.standO2./((profiledata.standSigma0 + 1000)/1000); %umol/L to umol/kg, normalized to density at surface (sigma0)
%Remove zero values (found in mission 1)
tol = 150; %remove values less than tolerance
[len, wid] = size(profiledata.standO2_raw);
for i = 1:len
    indbad = find(profiledata.standO2_raw(i,:) < tol);
    profiledata.standO2_raw(i,indbad) = NaN;
end
profiledata.standAOU_raw = profiledata.standO2equil - profiledata.standO2_raw; %apparent oxygen utilization (umol/kg)
profiledata.standO2biosat_raw = (profiledata.standO2_raw./profiledata.standO2equil - 1)*100; %O2 % biological supersaturation


%% Calibrate oxygen by prescribing that surface AOU match World Ocean Atlas climatology

%Calculate day-night index in order to only use nighttime data for AOU
%calibration
    UTCoffset = 0; %time is in UTC
    tol = 0.5; %number of hours before/after sunrise/sunset to count as daytime
[profiledata.dayind, profiledata.nightind] = indexDayNight(nanmean(profiledata.latlon_plot(:,1)), nanmean(profiledata.latlon_plot(:,2)), UTCoffset, profiledata.time, tol);

%Calculate offset to correct AOU to be zero at surface throughout
    %indnonan = find(isnan(profiledata.standAOU_raw(1, profiledata.nightind)) == 0);
    AOU_WOAclim_interp = interp1([WOAclim.dayinyr(1) - 30; WOAclim.dayinyr; WOAclim.dayinyr(end) + 30], [WOAclim.AOU(2,end) WOAclim.AOU(2,:) WOAclim.AOU(2,1)], profiledata.day(profiledata.nightind));
    offset = nanmean(AOU_WOAclim_interp - profiledata.standAOU_raw(1, profiledata.nightind)');
    residual = profiledata.standAOU_raw(1, profiledata.nightind)' + offset - AOU_WOAclim_interp;
    profiledata.residual_mean = nanmean(abs(residual));

if diagnostic_plots == 1
    figure(2); clf
    set(gcf,'color','w')
    x0=0; y0=0;
    width=20; height=17;
    set(gcf,'units','centimeters','position',[x0,y0,width,height])
        subplot(311)
    plot(profiledata.time, profiledata.standAOU_raw(1,:), 'k.'); hold on;
    plot(profiledata.time(profiledata.nightind), profiledata.standAOU_raw(1, profiledata.nightind), 'b.','markersize', 15); hold on;
    plot(profiledata.time(profiledata.nightind), AOU_WOAclim_interp, 'r--'); hold on;
    datetick('x'); title('Surface uncorr. AOU over mission')
    legend('All surface data','Nighttime surface data','World Ocean Atlas climatology','location','northwest')
        subplot(312)
    plot(profiledata.time, profiledata.standO2biosat_raw(1,:), 'k.'); hold on;
    plot(profiledata.time(profiledata.nightind), profiledata.standO2biosat_raw(1, profiledata.nightind), 'b.','markersize', 15); hold on;
    plot(profiledata.time(profiledata.nightind), zeros(length(profiledata.nightind)), 'k--'); hold on;
    datetick('x'); title('Surface uncorr. O_2 biosat over mission')
    legend('All surface data','Nighttime surface data','location','southwest')
        subplot(313)
    plot(profiledata.time(profiledata.nightind), residual, 'b.','markersize', 15); hold on;
    plot(profiledata.time(profiledata.nightind), zeros(length(profiledata.nightind)), 'k--'); hold on;
    datetick('x'); title('AOU residual after offset correction')
    
    print(['OxygenSensorCal_WOAClimatology_Mission' num2str(mission)],'-dpng')
end

%% Determine gain and offset corrections, and use to correct O2 data
%Currently prescribing no drift correction because no clear signal in
%extracting residual from WOA climatology
profiledata.gain = (nanmean(profiledata.standO2_raw(:,1)) - offset)/nanmean(profiledata.standO2_raw(:,1)); %based on mean surface O2 conc and offset correction for AOU
profiledata.drift = 0;

    [numdepths, ~] = size(profiledata.standO2_raw);
profiledata.standO2_corr = profiledata.standO2_raw.*profiledata.gain + profiledata.drift.*repmat((profiledata.time - min(profiledata.time)), 1, numdepths)'; %correct for gain and drift
profiledata.standAOU_corr = profiledata.standO2equil - profiledata.standO2_corr; %apparent oxygen utilization (umol/kg)
profiledata.standO2biosat_corr = (profiledata.standO2_corr./profiledata.standO2equil - 1)*100; %O2 % biological supersaturation

%% Show how O2 was corrected to surface value at equilibrium

if diagnostic_plots == 1
    figure(3); clf
        subplot(211)
    histogram(profiledata.standAOU_corr(1,:)); hold on;
    title('Histogram of corrected surface AOU over mission')
    xlim([-10 10])
        subplot(212)
    histogram(profiledata.standO2biosat_corr(1,:)); hold on;
    title('Histogram of corrected surface O_2 biosat over mission')
    xlim([-4 4])
end

%% Calculate DIC and NO3 from glider data using 2007 AR07E MLR

profiledata.DIC_wMLR07 = B_07(1,1) + B_07(2,1)*profiledata.standT + B_07(3,1)*profiledata.standsalt + B_07(4,1)*profiledata.standO2_corr + B_07(5,1)*profiledata.standAOU_corr;
profiledata.NO3_wMLR07 = B_07(1,2) + B_07(2,2)*profiledata.standT + B_07(3,2)*profiledata.standsalt + B_07(4,2)*profiledata.standO2_corr + B_07(5,2)*profiledata.standAOU_corr;

%% Plot corrected night-only data

if diagnostic_plots == 1
figure(4); clf
set(gcf,'color','w')
x0=0; y0=0;
width=35; height=20;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

        subplot(421)
    imagesc(profiledata.standConservT(:, profiledata.nightind)); colorbar;
    title('Conservative temperature');
    ylabel('dbar')
        subplot(423)
    imagesc(profiledata.standSA(:, profiledata.nightind)); colorbar;
    title('Absolute salinity');
    ylabel('dbar')
            subplot(425)
    imagesc(profiledata.NO3_wMLR07(:, profiledata.nightind)); colorbar;
    title('Nitrate from AR07E-07 MLR (\mumol/kg)');
    ylabel('dbar')
    xlabel('Nighttime profile number')

    if mission == 1
        subplot(427)
    imagesc(profiledata.standCHL(:, profiledata.nightind)); colorbar;
    title('Chlorophyll'); caxis([0 2])
    ylabel('dbar')
    end

        subplot(422)
    imagesc(profiledata.standO2_corr(:, profiledata.nightind)); colorbar;
    title('Corrected O_2 concentration (\mumol/kg)');
    ylabel('dbar')
        subplot(424)
    imagesc(profiledata.standAOU_corr(:, profiledata.nightind)); colorbar;
    ylabel('dbar')
    title('Corrected AOU (\mumol/kg)');
        subplot(426)
    imagesc(profiledata.standO2biosat_corr(:, profiledata.nightind)); colorbar;
    title('Corrected O_2 biosat (%)');
    ylabel('dbar')
        subplot(428)
    imagesc(profiledata.DIC_wMLR07(:, profiledata.nightind)); colorbar;
    ylabel('dbar')
    title('DIC from AR07E-07 MLR (\mumol/kg)');
    xlabel('Nighttime profile number')

    print(['GliderSectionsBGC_Mission' num2str(mission)],'-dpng')
end



%% Save data output from each mission

if mission == 1
    mission1 = profiledata; clear profiledata;
elseif mission == 2
    mission2 = profiledata; clear profiledata;
elseif mission == 4
     mission4 = profiledata; clear profiledata;
elseif mission == 5
    mission5 = profiledata; clear profiledata;
elseif mission == 6
    mission6 = profiledata; clear profiledata;
end

end

close all


%% Plot location of glider mission over map of AR07E occupations
C = [nicecolor('rrrk'); nicecolor('ry'); nicecolor('gbykw'); nicecolor('bbc'); nicecolor('brm')];

figure(1); clf
    latminplot = 48; latmaxplot = 68;
    lonminplot = -55; lonmaxplot = 0;
    M = 15;
m_proj('lambert','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
m_plot(data_05.lon, data_05.lat,'b.','markersize',M); hold on; %AR07E 2005
m_plot(data_07.lon, data_07.lat,'c.','markersize',M-1); hold on; %AR07E 2007
m_plot(data_14.lon(ind), data_14.lat(ind),'r.','markersize',M-2); hold on; %AR07E 2014
m_plot(mission1.latlon_plot(:,2), mission1.latlon_plot(:,1),'.','color',nicecolor('k'),'markersize',M/2); hold on;
m_plot(mission2.latlon_plot(:,2), mission2.latlon_plot(:,1),'.','color',nicecolor('k'),'markersize',M/2); hold on;
m_plot(mission4.latlon_plot(:,2), mission4.latlon_plot(:,1),'.','color',nicecolor('k'),'markersize',M/2); hold on;
m_plot(mission5.latlon_plot(:,2), mission5.latlon_plot(:,1),'.','color',nicecolor('k'),'markersize',M/2); hold on;
m_plot(mission6.latlon_plot(:,2), mission6.latlon_plot(:,1),'.','color',nicecolor('k'),'markersize',M/2); hold on;
m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
legend('','AR07E 2005','AR07E 2007','AR07E 2014','WHOUC glider missions','location','southwest')


%% Illustrate O2 calibration across all missions

C = [nicecolor('rrrk'); nicecolor('ry'); nicecolor('gbykw'); nicecolor('bbc'); nicecolor('brm')];
M = 10;

figure(2); clf
    subplot(211)
plot(mission1.day(mission1.nightind), mission1.standAOU_raw(1, mission1.nightind), '.', 'color', C(1,:), 'markersize', M); hold on;
plot(mission2.day(mission2.nightind), mission2.standAOU_raw(1, mission2.nightind),  '.', 'color', C(2,:), 'markersize', M); hold on;
plot(mission4.day(mission4.nightind), mission4.standAOU_raw(1, mission4.nightind),  '.', 'color', C(3,:), 'markersize', M); hold on;
plot(mission5.day(mission5.nightind), mission5.standAOU_raw(1, mission5.nightind),  '.', 'color', C(4,:), 'markersize', M); hold on;
plot(mission6.day(mission6.nightind), mission6.standAOU_raw(1, mission6.nightind),  '.', 'color', C(5,:), 'markersize', M); hold on;
plot(WOAclim.dayinyr, WOAclim.AOU(2,:), 'k.', 'markersize', M*2.5); hold on;
xlabel('Julian day'); ylabel('AOU, \mumol/kg')
title('Raw AOU from glider surface measurements')
    subplot(212)
plot(mission1.day(mission1.nightind), mission1.standAOU_corr(1, mission1.nightind), '.', 'color', C(1,:), 'markersize', M); hold on;
plot(mission2.day(mission2.nightind), mission2.standAOU_corr(1, mission2.nightind),  '.', 'color', C(2,:), 'markersize', M); hold on;
plot(mission4.day(mission4.nightind), mission4.standAOU_corr(1, mission4.nightind),  '.', 'color', C(3,:), 'markersize', M); hold on;
plot(mission5.day(mission5.nightind), mission5.standAOU_corr(1, mission5.nightind),  '.', 'color', C(4,:), 'markersize', M); hold on;
plot(mission6.day(mission6.nightind), mission6.standAOU_corr(1, mission6.nightind),  '.', 'color', C(5,:), 'markersize', M); hold on;
plot(WOAclim.dayinyr, WOAclim.AOU(2,:), 'k.', 'markersize', M*2.5); hold on;
legend('Mission 1', 'Mission 2', 'Mission 4', 'Mission 5', 'Mission 6', 'World Ocean Atlas 2013 climatology','location','southwest')
xlabel('Julian day'); ylabel('AOU, \mumol/kg')
title('Corrected AOU from glider surface measurements')

