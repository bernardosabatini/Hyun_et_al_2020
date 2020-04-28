% this runs the ROC on only the collected data points (as opposed to on the
% bootstrapped distrubitions)catOne=wPop4DCNO;
catZero=wPop4DSaline;

mergeF=[catZero' catOne'];
mergeL=[zeros(1, length(catZero)) ones(1, length(catOne))]';
[mergeX, mergeI]=sort(mergeF);
mergeL=mergeL(mergeI)';

np=length(mergeX);
outTPR=zeros(np, 1);
outFPR=zeros(np, 1);

nOne=length(catOne);
nZero=length(catZero);

for counter=1:np
    xNow=mergeX(counter);

    nOneHi=length(find(mergeL(counter:end)==1));
    nZeroHi=length(find(mergeL(counter:end)==0));

    %TP/(TP+FN) = TP/P
    TPR=nOneHi/nOne;
    %FP/(FN+TN) = FP/N
    FPR=nZeroHi/nZero;
    
    outTPR(counter)=TPR;
    outFPR(counter)=FPR;
    
end

figure
plot(outFPR, outTPR, 'LineWidth', 1, 'color', 'black')
set(gca, 'FontSize', 14)
title('4D CNO vs all Control')
    
-trapz(outFPR, outTPR)
    
    
    