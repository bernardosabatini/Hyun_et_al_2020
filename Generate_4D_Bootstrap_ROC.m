%% do the bootstrap and calculate confidence intervals
% there are 3 groups
%   1. CNO treated DREADD expressing animals
%   2. Drug vehicle treated DREADD expressing animals.  These are labeled
%       'saline' here even if the vehicle was something else
%   3. Control which is all the control data points including the vehicle
%       treated DREADD expressing, and both the CNO and vehicle treated
%       mCherry-only expressing animals

bt4DControl=bootstrp(10000, @mean, wPop4DControl);
[fi4DControl,xi4DControl] = ksdensity(bt4DControl);
ci4DControl=bootci(10000, @mean, wPop4DControl);

bt4DCNO=bootstrp(10000, @mean, wPop4DCNO);
[fi4DCNO,xi4DCNO] = ksdensity(bt4DCNO);
ci4DCNO=bootci(10000, @mean, wPop4DCNO);

bt4DSaline=bootstrp(10000, @mean, wPop4DSaline);
[fi4DSaline,xi4DSaline] = ksdensity(bt4DSaline);
ci4DSaline=bootci(10000, @mean, wPop4DSaline);

%% plot the results
figure
hold on

figMax=0.2;
ciY=figMax/4;

plot(xi4DControl, fi4DControl, 'color', 'black', 'LineWidth', 1)
plot(mean(wPop4DControl), 1.1*ciY, 'marker', 'X', 'LineStyle', 'none', 'MarkerEdgeColor', 'black', 'LineWidth', 1, 'MarkerSize', 8)  
errorbar(mean(wPop4DControl), 1.1*ciY, nan, nan, ci4DControl(1)-mean(wPop4DControl), ci4DControl(2)-mean(wPop4DControl), ...
    'color', 'black', 'LineWidth', 1)

plot(xi4DCNO, fi4DCNO, 'color', 'red', 'LineWidth', 1)
plot(mean(wPop4DCNO), ciY, 'marker', 'X', 'LineStyle', 'none', 'MarkerEdgeColor', 'red', 'LineWidth', 1, 'MarkerSize', 8)  
errorbar(mean(wPop4DCNO), ciY, nan, nan, ci4DCNO(1)-mean(wPop4DCNO), ci4DCNO(2)-mean(wPop4DCNO), ...
    'color', 'red', 'LineWidth', 1)

plot(xi4DSaline, fi4DSaline, 'color', 'green', 'LineWidth', 1)
plot(mean(wPop4DSaline), ciY, 'marker', 'X', 'LineStyle', 'none', 'MarkerEdgeColor', 'green', 'LineWidth', 1, 'MarkerSize', 8)  
errorbar(mean(wPop4DSaline), ciY, nan, nan, ci4DSaline(1)-mean(wPop4DSaline), ci4DSaline(2)-mean(wPop4DSaline), ...
    'color', 'green', 'LineWidth', 1)

set(gca, 'YLim', [0 figMax], 'XLim', [0 200])
set(gca, 'FontSize', 14)

%% Set up for the ROC anlaysis
catOneX=xi4DCNO;        % what is the first group
catOneF=fi4DCNO;        % this one gets the category label 1

catZeroX=xi4DControl;    % what is the other group
catZeroF=fi4DControl;    % this one gets the category label 2

dOne=catOneX(2)-catOneX(1);
dZero=catZeroX(2)-catZeroX(1);

mergeX=[catZeroX catOneX]';
mergeF=[dZero*catZeroF dOne*catOneF]';
mergeL=[zeros(1, length(catZeroX)) ones(1, length(catOneF))]';

xSort=sort(unique(mergeX));

np=length(xSort);
outTPR=zeros(np, 1);
outFPR=zeros(np, 1);

%% sweep through the thresholds
for counter=1:np
    xNow=xSort(counter);
    xUpToNow=find(mergeX<xNow);
    xMoreThanNow=find(mergeX>=xNow);
    
    if isempty(xUpToNow)
        nZeroLow=0;
        nOneLow=0;
    else
        nZeroLow=sum(mergeF((mergeL==0) & mergeX<xNow));
        nOneLow=sum(mergeF((mergeL==1) & mergeX<xNow));
    end        
    
    if isempty(xMoreThanNow)
        nZeroHi=0;
        nOneHi=0;
    else
        nZeroHi=sum(mergeF((mergeL==0) & mergeX>=xNow));
        nOneHi=sum(mergeF((mergeL==1) & mergeX>=xNow));
    end        
       
    %TP/(TP+FN) = TP/P
    TPR=nOneHi/(nOneHi+nOneLow);
    %FP/(FN+TN) = FP/N
    FPR=nZeroHi/(nZeroLow+nZeroHi);
    
    outTPR(counter)=TPR;
    outFPR(counter)=FPR;
    
end

%% plot the data
figure
plot(outFPR, outTPR, 'LineWidth', 1, 'color', 'black')
set(gca, 'FontSize', 14)
title('4D CNO vs all Control')
    
disp(['AUC is : ' num2str(-trapz(outFPR, outTPR))])
    
    
    