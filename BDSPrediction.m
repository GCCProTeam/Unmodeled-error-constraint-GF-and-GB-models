
% BDS prediction
m = 2;
count = 1;
for n = 4   
    load('GFBDS.mat');
    % Prediction using POLY
    for i = 1:14
    %   Input     
        str = strcat('BDS=C', num2str(i)); 
        eval(str);
    %   initialization  
        start = BDS(1,2);
        BDS(:,2) = (BDS(:,2)-start)./30 + 1;
    %   m:degree，n:window length    
        [YuCe,RMS] = PFPredict(BDS,m,n);
        YuCe(:, 2) = (YuCe(:, 2) - 1).*30 + start;
        str = strcat('CPredict', num2str(i), '=YuCe');
        eval(str); 
        BDSPredictedRMS(i,1) = RMS; 
        clear BDS;  
    end
    BDSPredictedRMS(BDSPredictedRMS(:, 1)==0,:)=[];
    Result(count,1) = n;
    Result(count,2) = mean(BDSPredictedRMS(:,1));
    count = count+1;
%     clear;
end
    