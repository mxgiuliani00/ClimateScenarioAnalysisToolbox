
% script for the analysis of climate change projections from EURO-CORDEX
% project on the Muzza irrigation district.

clear
clc

cd('/Users/matteo/Poli/didattica/2015-16_AESA/lab/lab01_IPCC')
where = pwd;

%% load and process netcdf files of precipitation and temperature over control period (1.1.2001-31.12.2005)
PRfiles = dir('pr*historical*.nc') ;
TASfiles = dir('tas*historical*.nc') ;

for i=1:length(PRfiles)
    disp(PRfiles(i).name)
    pr(i) = extractNetCDFdata( PRfiles(i).name, 'pr') ;
    disp(TASfiles(i).name)
    tas(i) = extractNetCDFdata( TASfiles(i).name, 'tas') ;
end
cd(where)

% assign structures
prH_2001_2005 = pr(1);
tasH_2001_2005 = tas(1);


%% visualizing time-snapshots: 9.12.2001
figure; plotNetCDFsnapshot(tasH_2001_2005, [2001,1,1], [2005,12,31], [2001,12,9], 'tas');
figure; plotNetCDFsnapshot(prH_2001_2005, [2001,1,1], [2005,12,31], [2001,12,9], 'pr');

% overlapping Italy shapefile 
cd ./shape
SHP = shaperead('ITA_adm0.shp');
cd (where)

% temperature
figure; plotNetCDFsnapshot(tasH_2001_2005, [2001,1,1], [2005,12,31], [2001,12,9], 'tas');
hold on; plot(SHP.X, SHP.Y, 'k');

% precipitation
figure; plotNetCDFsnapshot(prH_2001_2005, [2001,1,1], [2005,12,31], [2001,12,9], 'pr');
hold on; plot(SHP.X, SHP.Y, 'k');

%% RCP8.5 projection (1.1.2046-31.12-2050)-(1.1.2096-31.12.2100)
PRfiles = dir('pr*rcp85*.nc') ;
TASfiles = dir('tas*rcp85*.nc') ;

for i=1:length(PRfiles)
    disp(PRfiles(i).name)
    pr(i) = extractNetCDFdata( PRfiles(i).name, 'pr') ;
    disp(TASfiles(i).name)
    tas(i) = extractNetCDFdata( TASfiles(i).name, 'tas') ;
end
cd(where)

% assign structures
pr85_2046_2050 = pr(1);
pr85_2096_2100 = pr(2);
tas85_2046_2050 = tas(1);
tas85_2096_2100 = tas(2);

%% visualizing time-snapshots: 9.12.2001 vs 9.12.2050 vs 9.12.2100

% temperature
figure; 
subplot(131); plotNetCDFsnapshot(tasH_2001_2005, [2001,1,1], [2005,12,31], [2001,12,9], 'tas');
hold on; plot(SHP.X, SHP.Y, 'k');
subplot(132); plotNetCDFsnapshot(tas85_2046_2050, [2046,1,1], [2050,12,31], [2050,12,9], 'tas');
hold on; plot(SHP.X, SHP.Y, 'k');
subplot(133); plotNetCDFsnapshot(tas85_2096_2100, [2096,1,1], [2100,12,31], [2100,12,9], 'tas');
hold on; plot(SHP.X, SHP.Y, 'k');


% fix colormaps with the same limits 
x1 = plotNetCDFsnapshot(tasH_2001_2005, [2001,1,1], [2005,12,31], [2001,12,9], 'tas');
x2 = plotNetCDFsnapshot(tas85_2046_2050, [2046,1,1], [2050,12,31], [2050,12,9], 'tas');
x3 = plotNetCDFsnapshot(tas85_2096_2100, [2096,1,1], [2100,12,31], [2100,12,9], 'tas');
TT = [x1;x2;x3];
minT = min(TT(:));
maxT = max(TT(:));

figure; 
subplot(131); plotNetCDFsnapshot(tasH_2001_2005, [2001,1,1], [2005,12,31], [2001,12,9], 'tas');
hold on; plot(SHP.X, SHP.Y, 'k'); set(gca,'clim',[minT,maxT]); title('9/12/2001');
subplot(132); plotNetCDFsnapshot(tas85_2046_2050, [2046,1,1], [2050,12,31], [2050,12,9], 'tas');
hold on; plot(SHP.X, SHP.Y, 'k'); set(gca,'clim',[minT,maxT]); title('9/12/2050');
subplot(133); plotNetCDFsnapshot(tas85_2096_2100, [2096,1,1], [2100,12,31], [2100,12,9], 'tas');
hold on; plot(SHP.X, SHP.Y, 'k'); set(gca,'clim',[minT,maxT]); title('9/12/2100');


%% selection of location
lonMuzza = [9.3, 9.9];
latMuzza = [45.05, 45.5];

tasH_2001_2005_Muzza = cutNetCDFdata( tasH_2001_2005, lonMuzza, latMuzza );
tas85_2046_2050_Muzza = cutNetCDFdata( tas85_2046_2050, lonMuzza, latMuzza );
tas85_2096_2100_Muzza = cutNetCDFdata( tas85_2096_2100, lonMuzza, latMuzza );

prH_2001_2005_Muzza = cutNetCDFdata( prH_2001_2005, lonMuzza, latMuzza );
pr85_2046_2050_Muzza = cutNetCDFdata( pr85_2046_2050, lonMuzza, latMuzza );
pr85_2096_2100_Muzza = cutNetCDFdata( pr85_2096_2100, lonMuzza, latMuzza );
            
%% time-series analysis

% boxplot: central is median, the edges are 25th and 75th percentiles, the
% whiskers extend to extremes, outliers are plotted individually.
% outliers if they are larger (smaller) than Q3+1.5*(Q3-Q1), which
% corresponds to approximatively 2.7*sigma nad 99.3 coverage if the data
% are normally distributed
figure; 
subplot(121);boxplot([tasH_2001_2005_Muzza.value tas85_2046_2050_Muzza.value tas85_2096_2100_Muzza.value])
title('temperature'); set(gca,'XTickLabel', {'2001-2005', '2046-2050', '2096-2100'}); ylabel('°C')
subplot(122);boxplot([prH_2001_2005_Muzza.value pr85_2046_2050_Muzza.value pr85_2096_2100_Muzza.value])
title('precipitation'); set(gca,'XTickLabel', {'2001-2005', '2046-2050', '2096-2100'}); ylabel('mm/d')

% cyclostationary average
tasH_m = moving_average( tasH_2001_2005_Muzza.value, 365, 15 ) ; 
tas85_m46 = moving_average( tas85_2046_2050_Muzza.value, 365, 15 ) ;
tas85_m96 = moving_average( tas85_2096_2100_Muzza.value, 365, 15 ) ;
prH_m = moving_average( prH_2001_2005_Muzza.value, 365, 15 ) ; 
pr85_m46 = moving_average( pr85_2046_2050_Muzza.value, 365, 15 ) ;
pr85_m96 = moving_average( pr85_2096_2100_Muzza.value, 365, 15 ) ;

figure; 
subplot(121); plot( [tasH_m, tas85_m46, tas85_m96] ) ;
axis([1 365 0 30]); ylabel('temperature (°C)'); 
legend('hist', '2046-2050', '2096-2100');
subplot(122); plot( [prH_m, pr85_m46, pr85_m96] ) ;
axis([1 365 0 4]); ylabel('precipitation (mm/d)'); 
legend('hist', '2046-2050', '2096-2100');

% trend analysis
TT = [reshape(tasH_2001_2005_Muzza.value,365,5),...
    reshape(tas85_2046_2050_Muzza.value,365,5),...
    reshape(tas85_2096_2100_Muzza.value,365,5) ] ;

PP = [reshape(prH_2001_2005_Muzza.value,365,5),...
    reshape(pr85_2046_2050_Muzza.value,365,5),...
    reshape(pr85_2096_2100_Muzza.value,365,5) ] ;

[mT,hT]= MASH_plot(TT, 2001, 7, 2, 'temp');
[mP,hP]= MASH_plot(PP, 2001, 7, 3, 'prec');
