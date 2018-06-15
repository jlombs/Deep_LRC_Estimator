% Timing test comparison between LRCEstimator and hivePts

mu = 7*[40,30,20,10];
nu = mu;
lambda = 7*[65,55,46,34];

LRCTrue = 133792;

eps1 = .99;
eps2 = 1;
epsStore = [eps1,eps2];

numTrials = 60;

f1 = @() LRCEstimator(mu,nu,lambda,eps1,false);
f2 = @() LRCEstimatorC(mu,nu,lambda,eps2,false);
fhandles = [{f1},{f2}];

timeStore = [0,0];
LRCStores = zeros(2,numTrials);

%%{
fprintf(' Timings \n')
parfor i = 1:2
    
    timeStore(i) = timeit(fhandles{i});
    
end
%}
%timeStore(:) = 1;

fprintf('Accuracy \n')
parfor i = 1:numTrials

    LRCStores(1,i) =  LRCEstimator(mu,nu,lambda,eps1,false);
    
end
fprintf('Finished First\n')

parfor i = 1:numTrials
        
    LRCStores(2,i) =  LRCEstimatorC(mu,nu,lambda,eps2,false);
    
end

accuracyStore = abs(mean(LRCStores,2)-LRCTrue)/LRCTrue;
stdStore = (std(LRCStores,[],2)/sqrt(numTrials))./mean(LRCStores,2);

for i = 1:2
    fprintf('Algorithm %d with Eps = %2.1f Achieved %3.2f%% Accuracy with %3.2f%% fractional error in an Average of %4.3f Seconds \n',i,epsStore(i),100*accuracyStore(i),100*stdStore(i),...
        timeStore(i))
end


%{

Algorithm 1 with Eps = 0.4 Achieved 2.72% Error in an Average of 2.138 Seconds 
Algorithm 2 with Eps = 0.2 Achieved 5.98% Error in an Average of 2.006 Seconds 

Algorithm 1 with Eps = 0.3 Achieved .82% Error in an Average of 3.883 Seconds 
Algorithm 2 with Eps = 0.3 Achieved 2.93% Error in an Average of 1.348 Seconds 

Algorithm 1 with Eps = 0.5 Achieved 3.30% Error in an Average of 1.709 Seconds 
Algorithm 2 with Eps = 0.4 Achieved 1.18% Error in an Average of 1.149 Seconds 

Algorithm 1 with Eps = 0.9 Achieved 12.19% Error in an Average of 0.880 Seconds 
Algorithm 2 with Eps = 0.9 Achieved 0.85% Error in an Average of 0.933 Seconds 

Algorithm 2 with Eps = 0.1 Achieved 3.22% Error in an Average of 4.672 Seconds 

%}