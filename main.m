%main.m - Makes all data generation calls for desired days, prns parameters
%and contains some visualization methods in both 2D and 3D of angle and DCB
%variation correlation. Rudimentary data analysis includes binning data
%into (magnitude-indicated) angular increments and rotational slices in the angular x, y antenna space 
%
% Inputs (start of program)
%    doy - desired day of SV history
%    year - year of desired SV history
%    prns - array of desired SVs
%
% Outputs:
%    Locally saved SV orbit data (svOrbitClock), SV angle and DCB variation
%    data (generate_data) and selected 2D/3D plots
%
% Other files required: navsu repo, svOrbitClock.m, load_sv_posvelr.m,
% generate_data
%
% See also: svOrbitClock,  load_sv_posvelr, generate_data

% Author: Martin Freeman

clc;
clear;

cd('C:\Users\martin\Documents\MATLAB\DCB');
% misc PRN lists %
NO_CA_PRNS = [27, 9, 32, 26, 10, 30];
% IIRM_PRNS = [5,7,12,15,17,29,31];
% IIRM_SVNS = [50,48,58,55,53,57,52];
% JPL_PRNS = [17, 31, 12, 15, 29, 7];
% O_PRNS = [17, 31, 12, 29, 7];
% day_set1 = [32,60,91,121,152,182,213,244,274,305,335];

days = 2:363;
year = 2018;
prns = NO_CA_PRNS;

% generate a svn list
IIRM_SVNS = [];
for i=1:length(IIRM_PRNS)
    IIRM_SVNS = [IIRM_SVNS, navsu.svprn.prn2svn(IIRM_PRNS(i),315964800)];
end

%% process data
for i=1:length(days)
    generate_data(year, days(i), prns)
    fprintf('Day %d angle, dcb difference saved\n', days(i));
end

%% plot every PRN for every day
% cd('C:\Users\martin\Documents\MATLAB\DCB\SV_ANGLES_VAR\');
% for i=1:length(days)
%     for j=1:length(prns)
%         
%         doy = days(i);
%         prn = prns(j);
%         subplot(length(prns),length(days),length(days)*(j-1)+i)
%         
%         fdir = [num2str(year),'_',num2str(doy),'\',num2str(prn),'\SV_MTX.mat'];
%         SV_MTX = load(fdir);
%         SV_MTX = SV_MTX.SV_MTX;
%         plot3(SV_MTX(:,6),SV_MTX(:,7),SV_MTX(:,5),'b.');
%         title(['DCB Variation for PRN ',num2str(prns(j))]);
%         xlabel('Azimuth (x axis), deg');
%         ylabel('Elevation (y axis), deg');
%         zlabel('DCB Variation');
%         grid on;
%         hold on;
%         plot3([0,0],[0,0],[-7,7],'Color','k','linewidth',1);
%         view(90,0)
%     end
% end



%% grab aggragate, individual (all days, all SVs)
cd('C:\Users\martin\Documents\MATLAB\DCB\SV_ANGLES_VAR\');
AGG_SV_MTX = [];
IND_SV_MTX = {};
for i=1:length(days)
    for j=1:length(prns)
        doy = days(i);
        prn = prns(j);        
        fdir = [num2str(year),'_',num2str(doy),'\',num2str(prn),'\SV_MTX.mat'];
        SV_MTX = load(fdir);
        SV_MTX = SV_MTX.SV_MTX;
        IND_SV_MTX{i,j} = [SV_MTX(:,6:7),-SV_MTX(:,5),SV_MTX(:,8)]; % grab individual satellite
        AGG_SV_MTX = [AGG_SV_MTX; [SV_MTX(:,6:7),-SV_MTX(:,5),SV_MTX(:,8)]];
    end
end

%% bin analysis of SVs, 2D

incvec = linspace(0.1,2,20); % bin definition by bounds and increment size
figure;
for m=1:length(incvec)
    % construct new angle, variation array for bin angle value incvec(m)
    BINNED_SV_MTX = {};
    BINNED_FILTERED_SV_MTX = {};
    for i=1:length(prns)
        SV_MTX = cat(1,IND_SV_MTX{:,i});

    %     %log spacing
    %     theta_z_bins = flip(([0,logspace(-1,1.15,50)]-14.25)*(-1));
    %     n_bins = length(theta_z_bins)-1;
    
        %linear spacing
        max_theta_z = max(AGG_SV_MTX(:,4));
        theta_z_inc = incvec(m);
        n_bins = (max_theta_z-0)/theta_z_inc;
        theta_z_bins = linspace(0,max_theta_z, n_bins+1);
        for j=1:n_bins
            idx = find(SV_MTX(:,4) >= theta_z_bins(j) & SV_MTX(:,4) < theta_z_bins(j+1));
            BINNED_SV_MTX{i,j} = SV_MTX(idx,:);
            TF = isoutlier(BINNED_SV_MTX{i,j}(:,3),'mean','ThresholdFactor',3);
            BINNED_FILTERED_SV_MTX{i,j} = BINNED_SV_MTX{i,j}(~TF,:);
        end

    end


    % add to plot of individual SVs, binned by angle 
    
    for i=1:length(prns)
    %     figure(1);
    %     subplot(1,1,i);
    %     SV_MTX = cat(1,IND_SV_MTX{:,i});
    %     plot3(SV_MTX(:,1),SV_MTX(:,2),SV_MTX(:,3),'r.');
    %     title(['Ind. SV Code-Phase Bias SVN ',num2str(IIRM_SVNS(i)), ' Unprocessed'])
    %     xlabel('Azimuth (x axis), deg');
    %     ylabel('Elevation (y axis), deg');
    %     zlabel('Code-Phase Residual');
    %     grid on;
    %     hold on;
    %     plot3([0,0],[0,0],[-7,7],'Color','k','linewidth',1);
    %     hold on;
    %     ylim([-15,15]);
    %     view(0,0);
    %     
    %     figure(2);
    %     subplot(1,1,i);
    %     SV_MTX = sortrows(SV_MTX);
    %     plot(SV_MTX(:,4),SV_MTX(:,3),'r.');
    %     title(['CP Res and Nadir Angle, SVN ',num2str(IIRM_SVNS(i))]);
    %     xlabel('Angle, deg');
    %     ylabel('CP Res, m');
    %     grid on;
    %     ylim([-15,15]);
    %     
    %     figure(3);
    %     subplot(1,1,i);
    %     window_size = 50;
    %     TF = isoutlier(SV_MTX(:,3),'movmean',window_size);
    %     SV_MTX_FILTERED = SV_MTX(~TF,:);
    %     plot(SV_MTX_FILTERED(:,4),SV_MTX_FILTERED(:,3),'b.');
    %     title(['Filtered CP Res and Nadir Angle, SVN ',num2str(IIRM_SVNS(i))]);    
    %     xlabel('Angle, deg');
    %     ylabel('CP Res, m');
    %     grid on;
    %     ylim([-15,15]);
    %     
    %     figure(4);
    %     subplot(1,1,i);
    %     for j=1:n_bins
    %         mtx =  BINNED_SV_MTX{i,j};
    %         th = mean([theta_z_bins(j+1),theta_z_bins(j)]);
    %         plot(th,mean(mtx(:,3)),'b.');
    %         hold on;
    %     end
    %     title(['Unfiltered, Avg Binned CP Res and Nadir Angle, SVN ',num2str(IIRM_SVNS(i))]);    
    %     xlabel('Angle, deg');
    %     ylabel('CP Res, m');
    %     grid on;
    % %     ylim([-15,15]);
    %     
      
        subplot(4,5,m);
        th = [];
        pt = [];
        deg = 3;
        for j=1:n_bins
            %mtx =  BINNED_SV_MTX{i,j}; % no rejection
            mtx =  BINNED_FILTERED_SV_MTX{i,j};
            th = [th, mean([theta_z_bins(j+1),theta_z_bins(j)])];
            pt = [pt, mean(mtx(:,3))];
        end
        
        %fit polynomial to data
        %p = polyfix(th,pt,deg,th(1),pt(1)); 
        p = polyfit(th,pt,deg);
        f1 = polyval(p,th);
        SVN_SINGLE = navsu.svprn.prn2svn(prns(i),315964800);
        plot(th,pt,'b.');
        hold on;
        plot(th,f1,'g-','LineWidth',2);
        title({'Mean Res and Elevation',['SVN ',num2str(SVN_SINGLE),', ',num2str(length(days)),' Days, ',num2str(incvec(m)),' deg bin']});    
        xlabel('Angle, deg');
        ylabel('CP Res, m');
        grid on;
        legend('mean value',['deg ',num2str(deg),' poly fit']);
        xlim([0,15]);
        
    %     ylim([-15,15]);
    %     
    end
end


%% 3D Plots

% plot aggraggate data, unprocessed
figure;
plot3(AGG_SV_MTX(:,1),AGG_SV_MTX(:,2),AGG_SV_MTX(:,3),'r.');
title('Aggragate DCB Variation, unprocessed');
xlabel('Azimuth (x axis), deg');
ylabel('Elevation (y axis), deg');
zlabel('DCB Variation');
grid on;
hold on;
plot3([0,0],[0,0],[-7,7],'Color','k','linewidth',1);
view(0,0);


% generate and plot sliced, binned data in 3D. n slices around antenna and bin by theta_z

%sort into slices
[theta,rho] = cart2pol(AGG_SV_MTX(:,1),AGG_SV_MTX(:,2)); theta = wrapTo360(rad2deg(theta));
n_slices = 16;
ang_slices = linspace(0,360,n_slices+1);
SLICED_SV_MTX = {};
for i=1:n_slices
    idx = find(theta > ang_slices(i) & theta < ang_slices(i+1));
    SLICED_SV_MTX{i} = AGG_SV_MTX(idx,:);
end

%sort theta_z into bins of size theta_z_inc
binning_type = 1; %linear
switch binning_type

    case 1
        % linear binning: sort into theta_z bins
        max_theta_z = max(AGG_SV_MTX(:,4)); % observed value, for simplicity
        theta_z_inc = 0.3; % bin size for theta_z
        n_bins = (max_theta_z-0)/theta_z_inc;
        theta_z_bins = linspace(0,max_theta_z, n_bins+1);

        BINNED_SV_MTX = {};
        for i=1:n_slices
            for j=1:n_bins
                idx = find(SLICED_SV_MTX{i}(:,4) > theta_z_bins(j) & SLICED_SV_MTX{i}(:,4) < theta_z_bins(j+1));
                BINNED_SV_MTX{i,j} = SLICED_SV_MTX{i}(idx,:);
            end
        end
end


%optional filter step

for i=1:n_slices
    for j=1:n_bins
        mtx =  BINNED_SV_MTX{i,j};
        TF = isoutlier(mtx(:,3),'mean');
        BINNED_SV_MTX{i,j} = mtx(~TF,:);
    end
end

%plot mean of bins by slice. midpoint of theta_z increment used
MEAN_SV_MTX_BINNED = {};
figure;
for i=1:n_slices
    subplot(4,4,i);
    for j=1:n_bins
        mtx =  BINNED_SV_MTX{i,j};
        th = mean([theta_z_bins(j+1),theta_z_bins(j)]);
        sl_th = mean([ang_slices(i+1),ang_slices(i)]);
        MEAN_SV_MTX_BINNED{i,j} = [sl_th, th, mean(mtx(:,3))];
        plot(th,mean(mtx(:,3)),'r.');
        hold on;
    end
    title(['mean DCB variation with \theta_z, Slice: ',num2str(i)]);
    xlabel('\theta_z');
    ylabel('mean DCB variation (m)');
    grid on;
    ylim([-0.2,0.2]);
end

%re-aggragate into binned, 3D map
AGG_SV_MTX_BINNED = vertcat(MEAN_SV_MTX_BINNED{:,:});
[x,y] = pol2cart(AGG_SV_MTX_BINNED(:,1),AGG_SV_MTX_BINNED(:,2));
AGG_SV_MTX_BINNED = [x,y,AGG_SV_MTX_BINNED(:,3)];


figure;
subplot(1,2,1);
plotAggragate(AGG_SV_MTX_BINNED);
title('Aggragate CP Res, binned filtered');
subplot(1,2,2);
plotAggragateCubic(AGG_SV_MTX_BINNED);
title('Aggragate CP Res, Cubic fit, binned filtered');



%% plot variation with LOS vector Z axis angle ( > 0)
figure;
subplot(1,2,1);
AGG_SV_MTX = sortrows(AGG_SV_MTX); %x ang, y ang, dcb var, z ang
plot(AGG_SV_MTX(:,4),AGG_SV_MTX(:,3),'r.');
title('DCB Variation and Angle with Nadir');
xlabel('Angle, deg');
ylabel('DCB Variation, m');
grid on;
ylim([-15,15]);

subplot(1,2,2);
window_size = 25;
TF = isoutlier(AGG_SV_MTX(:,3),'movmean',window_size);
AGG_SV_MTX_FILTERED = AGG_SV_MTX(~TF,:);
plot(AGG_SV_MTX_FILTERED(:,4),AGG_SV_MTX_FILTERED(:,3),'r.');
title('Filtered DCB Variation and Angle with Nadir');
xlabel('Angle, deg');
ylabel('DCB Variation, m');
grid on;
ylim([-15,15]);


% misc helper functions
function plotAggragate(AGG_SV_MTX)
    xv = linspace(min(AGG_SV_MTX(:,1)),max(AGG_SV_MTX(:,1)),500);
    yv = linspace(min(AGG_SV_MTX(:,2)),max(AGG_SV_MTX(:,2)),500);
    [X,Y] = meshgrid(xv,yv);
    Z = griddata(AGG_SV_MTX(:,1),AGG_SV_MTX(:,2),AGG_SV_MTX(:,3),X,Y);
    mesh(X, Y, Z);
    title('Aggragate DCB Variation');
    xlabel('Azimuth (x axis), deg');
    ylabel('Elevation (y axis), deg');
    zlabel('DCB Variation');
    grid on;
    hold on;
    plot3([0,0],[0,0],[-.5,.5],'Color','k','linewidth',1);
    view(0,0);
end
function plotAggragateCubic(AGG_SV_MTX)
    xv = linspace(min(AGG_SV_MTX(:,1)),max(AGG_SV_MTX(:,1)),500);
    yv = linspace(min(AGG_SV_MTX(:,2)),max(AGG_SV_MTX(:,2)),500);
    [X,Y] = meshgrid(xv,yv);
    Z = griddata(AGG_SV_MTX(:,1),AGG_SV_MTX(:,2),AGG_SV_MTX(:,3),X,Y,'cubic');
    mesh(X, Y, Z);
    title('Aggragate DCB Variation');
    xlabel('Azimuth (x axis), deg');
    ylabel('Elevation (y axis), deg');
    zlabel('DCB Variation');
    grid on;
    hold on;
    plot3([0,0],[0,0],[-.5,.5],'Color','k','linewidth',1);
    view(0,0);
end