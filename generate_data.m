function generate_data(year, doy, prns)
%generate_data(year, doy, prns) - Sorts kaz data into structured,
%simplified, convenient arrays. Pulls IGS station solutions, calls orbit
%data generation script (svOrbitClock) to combine all necessary SV
%information and perform calculations/transformations to obtain SV_MTX 
%output (PRN, STN POS, DCB VAR, ANGLES X Y)
%
% Syntax: generate_data(year, doy, prns)
%
% Inputs:
%    doy - desired day of SV history
%    year - year of desired SV history
%    prns - array of desired SVs
%
% Outputs:
%    none
%
% Other files required: navsu repo, svOrbitClock.m, load_sv_posvelr.m
%
% See also: svOrbitClock,  load_sv_posvelr, main

% Author: Martin Freeman

% fix day string for later access
doy = num2str(doy);
while length(doy) < 3
    doy = ['0',doy];
end

%% load the data

% load kaz data set for given day, year
% data format: time index, meas type code/carrier, elevation angle, PRN, constellation
% index, predicted measurement, L1 C/N0, L2 C/N0, station index, Signal pair
% index (L1P-L2P/L1CA-L5), residual
cd 'C:\Users\martin\Documents\MATLAB\DCB\kaz_data'
data_name = ['clkDcbOutFull_',num2str(year),doy,'.mat'];
load(data_name);
cd 'C:\Users\martin\Documents\MATLAB\DCB'

%% organize into arrays with only useful data elements

% rm carrier phase
if length(dataFull.measInfoSave) > length(dataFull.b2)
    dataFull.measInfoSave = dataFull.measInfoSave(1:length(dataFull.b2),:);
else
    dataFull.b2 = dataFull.b2(1:length(dataFull.measInfoSave),:);
end
measInfo = [dataFull.measInfoSave, dataFull.b2];
idx = any(measInfo(:,10) > 0, 2);
L1L2Info = measInfo(~idx,:);
L1L5Info = measInfo(idx,:);

%match time index, PRN, station index
L1L2MatchMtx = [L1L2Info(:,1),L1L2Info(:,4), L1L2Info(:,9)]; 
L1L5MatchMtx = [L1L5Info(:,1),L1L5Info(:,4), L1L5Info(:,9)]; 
resultIdxA = ismember(L1L5MatchMtx, L1L2MatchMtx,'rows');
L1L5Info(~resultIdxA,:) = [];
resultIdxB = ismember(L1L2MatchMtx, L1L5MatchMtx,'rows');
L1L2Info(~resultIdxB,:) = [];

%make first three columns time index, PRN, station index for convenience
temp = [L1L2Info(:,4),L1L2Info(:,9),L1L5Info(:,4),L1L5Info(:,9)];
L1L2Info(:,4) = L1L2Info(:,2); L1L2Info(:,9) = L1L2Info(:,3); 
L1L5Info(:,4) = L1L5Info(:,2); L1L5Info(:,9) = L1L5Info(:,3); 
L1L2Info(:,2:3) = temp(:,1:2); L1L5Info(:,2:3) = temp(:,3:4); 
L1L2Info = sortrows(L1L2Info); L1L5Info = sortrows(L1L5Info);
clear temp L1L2mtx L1L5mtx L1L2MatchMtx L1L5MatchMtx idx measInfo resultIdxA resultIdxB

%% get the IGS station solutions

downloadChoice = 4; % 4 = IGS Station Solutions
year = str2num(data_name(15:18));
doy = str2num(data_name(19:21));

% get adjusted woy
jdi = navsu.time.doy2jd(year,doy);
[gpsWeek,gpsTow] = navsu.time.jd2gps(jdi);
gpsDow = floor(gpsTow/86400);
% Initial week of year
[gpsWeek0,tow0] = navsu.time.jd2gps(navsu.time.cal2jd(year,1,1));
woy = gpsWeek-gpsWeek0+1;
woy = num2str(woy);
while length(woy) < 2
    woy = ['0',woy];
end

% pulling station solutions manually before ftp fix
foldername = ['C:\Users\martin\Documents\MATLAB\DCB\data\precise-daily\',data_name(15:18) ...
,'\',data_name(19:21),'\'];
% 
% cd_comm = ['mkdir C:\Users\martin\Documents\MATLAB\DCB\data\precise-daily\',data_name(15:18) ...
% ,'\',data_name(19:21),'\'];
% status = system(cd_comm);
% cd_comm = ['cd C:\Users\martin\Documents\MATLAB\DCB\data\precise-daily\',data_name(15:18) ...
% ,'\',data_name(19:21),''];
% status = system(cd_comm);
% dl_comm = ['wget --ftp-user anonymous --ftp-password emmerich@stanford.edu ftps://gdc.cddis.eosdis.nasa.gov/gps/products/', ... 
%     num2str(gpsWeek),'/IGS18P',woy,'_all.ssc.Z ', '-P ',foldername];
% 
% 
% status = system(dl_comm);
% 
% 
% 
cd(foldername)
filename = ['IGS18P',woy,'_all.ssc.Z'];
% 
% path_to_7z='"C:\Program Files\7-Zip\7z.exe"'; % adapt to your path
% str = [path_to_7z, ' e ', filename]; % note the blank
% system(str);


% pulling solutions with working ftp
navsu.ftp.download(downloadChoice,year,doy); %done manually atm
fileList = dir(['C:\Users\martin\Documents\MATLAB\DCB\data\precise-daily\',data_name(15:18) ...
,'\',data_name(19:21),'\*.ssc']);
filename = ['C:\Users\martin\Documents\MATLAB\DCB\data\precise-daily\',data_name(15:18) ...
,'\',data_name(19:21),'\',fileList.name];

sscData = navsu.readfiles.parseSscFile(['IGS18P',woy,'_all.ssc']);

cd 'C:\Users\martin\Documents\MATLAB\DCB'


% check all station codes used in kaz dataset are pulled
statCodes = dataFull.statCodes;
statPos = [];
for i=1:length(statCodes)
    resultIdx = find(ismember(sscData.name, statCodes{i}));
    if isempty(resultIdx)
        resultIdx = find(ismember(L1L2Info(:,3),i));
        if ~isempty(resultIdx)
            disp('Missing a station. Check if all IGS solutions were loaded');
            continue
        end
        continue
    end
    statPos(i,1:3) = sscData.data(resultIdx,1:3);
end

%populate "information" array with station positions
for i=1:length(L1L2Info)
    InfoPos(i,:) = statPos(L1L2Info(i,3),:);
end


%% generate and locally save the sv trajectory data (pos, vel, rotn mtx R)

NO_CA_PRNS = [27, 9, 32, 26, 10, 30];
dt = 60; 
n_epochs = 1440;

% generate the sv trajectory data for all prns on doy, year
svOrbitClock(doy,year,NO_CA_PRNS, dt, n_epochs)

% access and load locally saved sv trajectory data into arrays
SV_POS_VEL = [];
SV_R = [];
for j=1:length(prns)
    cd('C:\Users\martin\Documents\MATLAB\DCB\');
    [svPos,svVel,R] = load_sv_posvelr(prns(j),year,doy);
    SV_POS_VEL(:,:,j) = [svPos.svPos, svVel.svVel];
    SV_R(:,:,:,j) = R.R;
end

clear svPos svVel R 
cd('C:\Users\martin\Documents\MATLAB\DCB\');

%% generate SV angles and DCB difference values
DCB_VAR = L1L2Info(:,11) - L1L5Info(:,11); % difference of L1L2, L1L5 residuals

%consolidate all relevant sv info into one mtx for operations
for i=1:length(prns)
    SV_MTX = [];
    prn = prns(i);
    prnIdx = find(ismember(L1L2Info(:,2),prn)); %make sure they're matched..
    SV_MTX(:,:) = [L1L2Info(prnIdx,1), InfoPos(prnIdx,:), DCB_VAR(prnIdx)];
    ANG_MTX = [];
    for j=1:length(SV_MTX) % transform data into antenna frame
        epoch = SV_MTX(j,1);
        VLOS_ECF = (SV_MTX(j,2:4)-SV_POS_VEL(epoch,1:3,i))/norm(SV_MTX(j,2:4)-SV_POS_VEL(epoch,1:3)); 
        VLOS_RSV = SV_R(:,:,epoch,i)'*VLOS_ECF';
        Z_ANG = acosd(dot(VLOS_RSV,[0,0,1])/(norm(VLOS_RSV)*norm([0,0,1])));
        ANG_MTX(j,:) = [atan2d(VLOS_RSV(1),VLOS_RSV(3)), atan2d(VLOS_RSV(2),VLOS_RSV(3)), Z_ANG];
    end
    SV_MTX = [SV_MTX, ANG_MTX]; %PRN, STN POS, DCB VAR, ANGLES X Y
    
    % save arrays with PRN, STN POS, DCB VAR, ANGLES X, Y
    cd('C:\Users\martin\Documents\MATLAB\DCB\SV_ANGLES_VAR\');
    fdir = [num2str(year),'_',num2str(doy),'\',num2str(prn),'\'];
    mkdir(fdir)
    cd(fdir)
    filename = ['C:\Users\martin\Documents\MATLAB\DCB\SV_ANGLES_VAR\',fdir];
    save([filename,'SV_MTX.mat'],'SV_MTX');
    cd('C:\Users\martin\Documents\MATLAB\DCB\');
end

end
