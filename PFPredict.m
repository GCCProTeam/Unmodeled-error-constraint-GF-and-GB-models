
% Fitting and prediction using POLY
function [YuCe,RMS] = PFPredict(SatData,m,n)

    D = SatData;
    windowcount = n;
    row = D(end, 2);
    count  = 1;
    countt = 1;    
    interval1 = 0;
    interval2 = 0;
    startcount = 1;
    endcount = windowcount;  

    while (1)       
        y2 = D(startcount:endcount, :);       
        y2(:,2) = (y2(:,2)-y2(1,2))+1;      
%       
        interval1 = y2(1, 2);
        interval2 = y2(end, 2);
        interval = interval2 - interval1;
        
        if interval > n-1
            countt = countt + 1;
            if countt == 2
                count = count - 1;
                YuCe(count, :)=[];
            end           
            if D(endcount + 1,2) == row
                break;
            end         
            startcount = startcount + 1;  
            endcount   = endcount + 1;
            continue;
        end
          
        p = polyfit(y2(:,2),y2(:,4),m);        
        YuCeT = polyval(p,(y2(end,2)+1));
        % Prediction
        YuCe(count,4) = YuCeT;
        YuCe(count,1:3) = D(endcount + 1,1:3);
        count = count + 1;
        countt = 1;      
        
        if D(endcount + 1,2) == row
            break;
        end       
        startcount = startcount + 1;  
        endcount   = endcount + 1; 

    end  
    [~, ia1, ib1] = intersect(D(:, 2),YuCe(:, 2));
    A = D(ia1(:, 1), 4);
    B = YuCe(ib1(:, 1), 4);
    
    RMS1 = (A - B);
    RMS1 = RMS1.^2;
    RMS2 = sum(RMS1);
    [r,~]= size(A);
    RMS3 = sqrt(RMS2/(r));
    RMS  = RMS3; 
end
