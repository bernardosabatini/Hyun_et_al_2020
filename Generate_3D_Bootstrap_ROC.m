%% do the bootstrap and calculate confidence intervals
% there are 3 groups
%   1. CNO treated DREADD expressing animals
%   2. Drug vehicle treated DREADD expressing animals.  These are labeled
%       'saline' here even if the vehicle was something else
%   3. Control which is all the control data points including the vehicle
%       treated DREADD expressing, and both the CNO and vehicle treated
%       mCherry-only expressing animals

bt3DControl=bootstrp(10000, @mean, wPop3DControl);
[fi3DControl,xi3DControl] = ksdensity(bt3DControl);
ci3DControl=bootci(10000, @mean, wPop3DControl);

bt3DCNO=bootstrp(10000, @mean, wPop3DCNO);
[fi3DCNO,xi3DCNO] = ksdensity(bt3DCNO);
ci3DCNO=bootci(10000, @mean, wPop3DCNO);

bt3DSaline=bootstrp(10000, @mean, wPop3DSaline);
[fi3DSaline,xi3DSaline] = ksdensity(bt3DSaline);
ci3DSaline=bootci(10000, @mean, wPop3DSaline);

%% plot the results
figure
hold on

figMax=0.06;
ciY=figMax/4;

plot(xi3DControl, fi3DControl, 'color', 'black', 'LineWidth', 1)
plot(mean(wPop3DControl), 1.1*ciY, 'marker', 'X', 'LineStyle', 'none', 'MarkerEdgeColor', 'black', 'LineWidth', 1, 'MarkerSize', 8)  
errorbar(mean(wPop3DControl), 1.1*ciY, nan, nan, ci3DControl(1)-mean(wPop3DControl), ci3DControl(2)-mean(wPop3DControl), ...
    'color', 'black', 'LineWidth', 1)

plot(xi3DCNO, fi3DCNO, 'color', 'red', 'LineWidth', 1)
plot(mean(wPop3DCNO), ciY, 'marker', 'X', 'LineStyle', 'none', 'MarkerEdgeColor', 'red', 'LineWidth', 1, 'MarkerSize', 8)  
errorbar(mean(wPop3DCNO), ciY, nan, nan, ci3DCNO(1)-mean(wPop3DCNO), ci3DCNO(2)-mean(wPop3DCNO), ...
    'color', 'red', 'LineWidth', 1)

plot(xi3DSaline, fi3DSaline, 'color', 'green', 'LineWidth', 1)
plot(mean(wPop3DSaline), ciY, 'marker', 'X', 'LineStyle', 'none', 'MarkerEdgeColor', 'green', 'LineWidth', 1, 'MarkerSize', 8)  
errorbar(mean(wPop3DSaline), ciY, nan, nan, ci3DSaline(1)-mean(wPop3DSaline), ci3DSaline(2)-mean(wPop3DSaline), ...
    'color', 'green', 'LineWidth', 1)

set(gca, 'YLim', [0 figMax], 'XLim', [0 500])
set(gca, 'FontSize', 14)

%% Set up for the ROC anlaysis
catZeroX=xi3DCNO;        % what is the first group
catZeroF=fi3DCNO;        % this one gets the category label 1

catOneX=xi3DControl;    % what is the other group
catOneF=fi3DControl;    % this one gets the category label 2

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
title('3D CNO vs all Control')
    
disp(['AUC is : ' num2str(-trapz(outFPR, outTPR))])
    
    
    