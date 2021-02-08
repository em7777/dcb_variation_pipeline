 function [svPos,svVel,R] = load_sv_posvelr(prn, year, doy)
    cd('C:\Users\martin\Documents\MATLAB\DCB\SV_STATES\');
    dir = [num2str(year),'_',num2str(doy),'\',num2str(prn),'\'];
    filename = ['C:\Users\martin\Documents\MATLAB\DCB\SV_STATES\',dir];
    svPos = load([filename,'svPos.mat']);
    svVel = load([filename,'svVel.mat']);
    R = load([filename,'svR.mat']);
 end
 
