function [status] = svOrbitClock(days, year, prns, dt, n_epochs)
%svOrbitClock(days,year, prns) - Pull data products of given PRNs for all days (days) in year (year), locally saves position, velocity and
%rotation mtx R history for each SV. Ideally this is run only once to locally store your desired SV history days/prns.
%
% Syntax:  [status] = svOrbitClock(days,year, prns)
%
% Inputs:
%    days - array of desired SV history days
%    year - year of desired SV history
%    prns - array of desired SVs
%    dt - time interval in seconds
%    n_epochs - number of epochs 
% Outputs:
%    status - (optional) confirmation of completion
%
%
% Other files required: navsu repo
%
% See also: generate_data,  main

% Author: Martin Freeman

%%%%%%% misc. day/prn arrays: %%%%%%%
% first_of_month_days = [32,60,91,121,152,182,213,244,274,305,335];
% NO_CA_PRNS = [27, 9, 32, 26, 10, 30];
% IIRM_PRNS = [5,7,12,15,17,29,31];
% days = 2:363;
% year = 2018;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=days(1):days(end)
    for j=1:length(prns)
        % Input vector for navsu.svOrbitClock telling it what constellations to
        % load.  First element is GPS, then GLONASS, Galileo, BDS, SBAS
        % respectively.  
        constUse = [1 1 0 0 0];

        % A configuration file should be specified- this tells the software where
        % to place all of the downloaded products.  It can be modeled after the
        % default.ini file that already exists. 
        configFile = 'config.ini';

        % Username and password file for NASA data/products download. See: 
        % [1] https://cddis.nasa.gov/Data_and_Derived_Products/CDDIS_Archive_Access.html
        % [2] https://cddis.nasa.gov/Data_and_Derived_Products/CreateNetrcFile.html
        netrcFile = 'C:\Users\martin\Documents\cddislogin.netrc';

        % Initialize the GNSS correction data object
        datai = navsu.svOrbitClock('constUse',[1 1 0 0 0],'netrcFile',netrcFile);

        % Year and day of year of interest- this is just arbitrary for this
        % example
        doy  = days(i);

        %% Choose an example time frame and satellite
        % Choose a time from the middle of the day (need to be able to interpolate
        % using data from both sides)

        % Start time is built from day of year information 
        epochStart = navsu.time.jd2epochs(navsu.time.doy2jd(year,doy)+0.5);  
        epochs = (epochStart:dt:(epochStart+n_epochs*60-60))'; 

        % Choosing an arbitrary GPS or GLONASS satellite 
        prn = prns(j);
        constInd = 1;  % GPS = 1, GLONASS = 2, GAL = 3

        % the PRN, constInd, and epoch inputs must all be the same length to the
        % orbit interpolation function
        prnVec = prn*ones(size(epochs));
        constIndVec = constInd*ones(size(epochs));

        %% Initialize Orbit and clock corrections
        disp('Downloading and/or parsing orbit and clock corrections')
        % Download and/or parse the orbital data
        datai.initOrbitData([year year year],[doy-1 doy doy+1]);

        % Download and/or parse the clock data
        %datai.initClockData([year year year],[doy-1 doy doy+1]);

        % Interpolate the precise orbit and clock data
        disp('Interpolating orbit and clock corrections to the desired epochs')

        % Orbit data: Get interpolated SV position, velocity, rotation
        % mtx R and sun position for yaw-steering frame transform
        [svPos,svVel] = datai.propagate(prnVec,constIndVec,epochs);
        [R, sunpos] = navsu.geo.svLocalFrame(svPos, epochs);

        % Clock data
        %svBias = datai.clock(prnVec,constIndVec,epochs);

        % Store data: Store orbit data (pos, vel, R) locally by sv, day
        cd('C:\Users\martin\Documents\MATLAB\DCB\SV_STATES\');
        dir = [num2str(year),'_',num2str(doy),'\',num2str(prn),'\'];
        mkdir(dir)
        cd(dir)
        filename = ['C:\Users\martin\Documents\MATLAB\DCB\SV_STATES\',dir];
        save([filename,'svPos.mat'],'svPos');
        save([filename,'svVel.mat'],'svVel');
        save([filename,'svR.mat'],'R');
        cd('C:\Users\martin\Documents\MATLAB\DCB\');

    end

    fprintf('Day %d orbit data saved\n', days(i));

end

return 
    
end




















