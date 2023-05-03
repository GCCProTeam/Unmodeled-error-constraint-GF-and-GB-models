
% GPS prediction
m = 2;
count = 1;
for n = 4
    load('GFGPS.mat');
    % Prediction using POLY
    for i = 1:32    
    %   Input    
        str = strcat('GPS=G', num2str(i)); 
        eval(str);
    %   initialization 
        start = GPS(1,2);
        GPS(:,2) = (GPS(:,2)-start)./30 + 1;
    %   m:degree，n:window length  
        [YuCe,RMS] = PFPredict(GPS,m,n);
        YuCe(:, 2) = (YuCe(:, 2) - 1).*30 + start;
        str = strcat('GPredict', num2str(i), '=YuCe');
        eval(str); 
        GPSPredictedRMS(i,1) = RMS; 
        clear GPS;  
    end
    GPSPredictedRMS(GPSPredictedRMS(:, 1)==0,:)=[];
    Result(count,1) = n;
    Result(count,2) = mean(GPSPredictedRMS(:,1));
    count = count+1;
%     clear;
end  
