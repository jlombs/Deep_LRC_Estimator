% LRC Accuracy plots, Cousins

close all
clear 
clc

%{
mu = [2 1 0];
nu = mu;
lambda = [3 2 1];
%}

%%{
mu = [40,30,20,10];
nu = mu;
lambda = [65,55,46,34];

%}

fileID = fopen(sprintf('./LRCoef/outputData%d.txt',numel(mu)),'r');
A = fscanf(fileID,'%d');
fclose(fileID);

d = A(1);
n = A(2:2:end)';
c = A(3:2:end)';

epsRange = [.1,.9];

LRCDataC = zeros(2,numel(n),numel(epsRange));
LRCData = zeros(2,numel(n),numel(epsRange));

numSamps = 60;

counter = 1;
for i = n
    
    fprintf('%d/%d \n',i,max(n))
    tmpC = zeros(1,numSamps);
    tmp = zeros(1,numSamps);
    
    counterE = 1;
    
    for eps = epsRange
    %{
        parfor j = 1:numSamps

            tmpC(j) = LRCEstimatorC(i*mu,i*nu,i*lambda,eps,false);

        end
        %}
        parfor j = 1:numSamps

            tmp(j) = LRCEstimator(i*mu,i*nu,i*lambda,eps,false);

        end

        LRCData(1,counter,counterE) = mean(tmp);
        LRCData(2,counter,counterE) = 1.96*std(tmp)/sqrt(numSamps);
        %{
        LRCDataC(1,counter,counterE) = mean(tmpC);
        LRCDataC(2,counter,counterE) = 1.96*std(tmpC)/sqrt(numSamps);
        %}
        counterE = counterE + 1;
        
    end
    
    counter = counter + 1;
    
end

%{
legendStorage = cell(1,3);
colors = ['r','b'];
figure

for i = 1:2
    
    errorbar(n,LRCDataC(1,:,i),LRCDataC(2,:,i),sprintf('-%s',colors(i)))
    legendStorage{i} = sprintf('Cousins Estimate with eps=%3.2f',epsRange(i));
    hold on
    
end

legendStorage{end} = 'Exact';

scatter(n,c,'k')
title([sprintf('Cousins Estimated LRC on %d Samples Compared to Known Value for Weight Multiples: \n [',numSamps),...
    sprintf('%d ',lambda),sprintf(']; ['),sprintf('%d ',nu),sprintf(']; ['),...
    sprintf('%d ',mu),sprintf(']')])
xlabel('Weight Multiple')
ylabel('LRC')
legend(legendStorage,'Location','SE')
%}

legendStorage = cell(1,3);
colors = ['r','b'];
figure
for i = 1:2
    
    errorbar(n,LRCData(1,:,i),LRCData(2,:,i),sprintf('-%s',colors(i)))
    legendStorage{i} = sprintf('CHAR Estimate with eps=%3.2f',epsRange(i));
    hold on
    
end

legendStorage{end} = 'Exact';

scatter(n,c,'k')
title([sprintf('CHAR Estimated LRC on %d Samples Compared to Known Value for Weight Multiples: \n [',numSamps),...
    sprintf('%d ',lambda),sprintf(']; ['),sprintf('%d ',nu),sprintf(']; ['),...
    sprintf('%d ',mu),sprintf(']')])
xlabel('Weight Multiple')
ylabel('LRC')
legend(legendStorage,'Location','SE')
