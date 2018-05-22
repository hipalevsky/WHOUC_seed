%% Read in AR07E bottle data for 2005 and 2007 sections and calculate MLRs for DIC and NO3
% Predictor variables are T, S, O2, and AOU
[data_05, B_05, residmean_05, R2_05] = MLR_fromGOSHIPdata('C:/Users/Hilary/Dropbox/WHOI Postdoc/OUC seed proposal/AR07E-2014/64PE20050907.csv', 2); %note that 2nd input variable is fudge factor to deal with different arrangement of columns
[data_07, B_07, residmean_07, R2_07] = MLR_fromGOSHIPdata('C:/Users/Hilary/Dropbox/WHOI Postdoc/OUC seed proposal/AR07E-2014/64PE20070830.csv', 0);
[data_14, B_14, residmean_14, R2_14] = MLR_fromGOSHIPdata_netcdf('C:/Users/Hilary/Dropbox/WHOI Postdoc/OUC seed proposal/AR07E-2014/sam_jr302_all.nc');

%% Test applicability of MLR from one cruise to predict values on the other

%Calculate DIC and NO3 in 2007 using 2005 MLRs
DIC_07_wMLR05 = B_05(1,1) + B_05(2,1)*data_07.T + B_05(3,1)*data_07.S + B_05(4,1)*data_07.O2 + B_05(5,1)*data_07.AOU;
NO3_07_wMLR05 = B_05(1,2) + B_05(2,2)*data_07.T + B_05(3,2)*data_07.S + B_05(4,2)*data_07.O2 + B_05(5,2)*data_07.AOU;
    ind = find(data_07.pres > 30 & data_07.lat > 55 & data_07.lon > -42.5); %remove surface values and outlier stations - this is the same as in function (could change to an input variable)
    residmean_07_wMLR05 = [nanmean(abs(DIC_07_wMLR05(ind) - data_07.DIC(ind))), nanmean(abs(NO3_07_wMLR05(ind) - data_07.NO3(ind)))];

%Calculate DIC and NO3 in 2005 using 2007 MLRs 
DIC_05_wMLR07 = B_07(1,1) + B_07(2,1)*data_05.T + B_07(3,1)*data_05.S + B_07(4,1)*data_05.O2 + B_07(5,1)*data_05.AOU;
NO3_05_wMLR07 = B_07(1,2) + B_07(2,2)*data_05.T + B_07(3,2)*data_05.S + B_07(4,2)*data_05.O2 + B_07(5,2)*data_05.AOU;
    ind = find(data_05.pres > 30 & data_05.lat > 55 & data_05.lon > -42.5); %remove surface values and outlier stations - this is the same as in function (could change to an input variable)
    residmean_05_wMLR07 = [nanmean(abs(DIC_05_wMLR07(ind) - data_05.DIC(ind))), nanmean(abs(NO3_05_wMLR07(ind) - data_05.NO3(ind)))];

%Calculate DIC and NO3 in 2014 using 2007 MLRs
DIC_14_wMLR07 = B_07(1,1) + B_07(2,1)*data_14.T + B_07(3,1)*data_14.S + B_07(4,1)*data_14.O2 + B_07(5,1)*data_14.AOU;
NO3_14_wMLR07 = B_07(1,2) + B_07(2,2)*data_14.T + B_07(3,2)*data_14.S + B_07(4,2)*data_14.O2 + B_07(5,2)*data_14.AOU;
    ind = find(data_14.pres > 30 & data_14.lat > 55 & data_14.lon > -42.5); %remove surface values and outlier stations - this is the same as in function (could change to an input variable)
    residmean_14_wMLR07 = [nanmean(abs(DIC_14_wMLR07(ind) - data_14.DIC(ind))), nanmean(abs(NO3_14_wMLR07(ind) - data_14.NO3(ind)))];
    
%%% Able to adequately determine sub-ML DIC and NO3 using MLR from a
%%% different year (though note that a long time difference will start to
%%% introduce anthropogenic increase in DIC). Also could try using this in the
%%% Irminger Sea or elsewhere along section...

