% LRC Accuracy plots, Cousins

close all
clear 
clc

%%{
mu = [2 1 0];
nu = mu;
lambda = [3 2 1];
%}

%{
mu = [40,30,20,10];
nu = mu;
lambda = [65,55,46,34];

%}

%{
mu = [40,30,20,10];
nu = mu;
lambda = [65,55,45,35];
%}

fileID = fopen(sprintf('./LRCoef/outputData%d.txt',numel(mu)),'r');
A = fscanf(fileID,'%d');
fclose(fileID);

d = A(1);
n = A(2:2:end)';
c = A(3:2:end)';

n = n(1:(end-2));
c = c(1:(end-2));

epsRange = [.2];

LRCData = zeros(2,numel(n),numel(epsRange));

numSamps = 30;

counter = 1;
for i = n
    
    fprintf('%d/%d \n',i,max(n))
    tmp = zeros(1,numSamps);
    
    counterE = 1;
    
    for eps = epsRange
    
        parfor j = 1:numSamps

            tmp(j) = LRCEstimatorC(i*mu,i*nu,i*lambda,eps,false);

        end

        LRCData(1,counter,counterE) = mean(tmp);
        LRCData(2,counter,counterE) = 1.96*std(tmp)/sqrt(numSamps);
        
        counterE = counterE + 1;
        
    end
    
    counter = counter + 1;
    
end

figure
legendStorage = cell(1,numel(epsRange)+1);
for i = 1:numel(epsRange)
    
    errorbar(n,LRCData(1,:,i),LRCData(2,:,i))
    legendStorage{i} = sprintf('Estimate with eps=%3.2f',epsRange(i));
    hold on
    
end
legendStorage{end} = 'Exact';

scatter(n,c,'k')
title([sprintf('Estimated LRC on %d Samples Compared to Known Value for Weight Multiples: \n [',numSamps),...
    sprintf('%d ',lambda),sprintf(']; ['),sprintf('%d ',nu),sprintf(']; ['),...
    sprintf('%d ',mu),sprintf(']')])
xlabel('Weight Multiple')
ylabel('LRC')
legend(legendStorage,'Location','SE')

figure
legendStorage = cell(1,numel(epsRange));
for i = 1:numel(epsRange)
    
    scatter(n,100*abs(c-LRCData(1,:,i))./c)
    legendStorage{i} = sprintf('Estimate with eps=%3.2f',epsRange(i));
    hold on
    
end

title([sprintf('Estimated LRC Percent Error on %d Samples Compared to Known Value for Weight Multiples: \n [',numSamps),...
    sprintf('%d ',lambda),sprintf(']; ['),sprintf('%d ',nu),sprintf(']; ['),...
    sprintf('%d ',mu),sprintf(']')])
xlabel('Weight Multiple')
ylabel('LRC')
legend(legendStorage,'Location','SE')

