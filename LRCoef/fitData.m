%Script to fit LRC data with non-neg poly coefficients

close all

fileID = fopen('outputData4.txt','r');
A = fscanf(fileID,'%d');
fclose(fileID);

d = A(1);
n = A(2:2:end);
c = A(3:2:end);
k = nchoosek(d-1,2);

%% Build 'LP' for polyfit of order k
pts = numel(n);
A = ones(pts,k+1);

for i = 1:pts
    A(i,2:end) = arrayfun(@(z) n(i)^z,1:k);
end

%% Solve LP
[LP,resnorm] = lsqnonneg(A,c);
LP = LP';

numSamp = 10*pts;
xp = linspace(n(1),n(end),numSamp);
fp = ones(1,numSamp);
for i = 1:numSamp
    fp(i) = sum(LP.*arrayfun(@(z) xp(i)^z,0:k));
end

figure
scatter(n,c)
hold on
plot(xp,fp,'r--')
title(sprintf('Leading Power Coefficient = %3.2f\n Residual Norm = %3.1e',LP(end),resnorm))
xlabel('Multiples of Weight Vector')
ylabel('LRCs')
legend(sprintf('Dimension = %d',d))

%% Noise Injection
% Want to study the robustness of the leading coefficient when the data is
% perturbed by gaussian noise.

noiseLevel = .05;
cNoise = c+c.*(noiseLevel*(randn(pts,1)));
%{
figure
scatter(n,c,'k')
hold on
scatter(n,cNoise,'b')
%}

[LP,resnorm,residual] = lsqnonneg(A,cNoise);
LP = LP';

numSamp = 10*pts;
xp = linspace(n(1),n(end),numSamp);
fp = ones(1,numSamp);
for i = 1:numSamp
    fp(i) = sum(LP.*arrayfun(@(z) xp(i)^z,0:k));
end

figure
scatter(n,c,'k')
hold on
scatter(n,cNoise,'b')
plot(xp,fp,'r--')
title(sprintf('Leading Power Coefficient = %3.2f\n Residual Norm = %3.1e',LP(end),resnorm))
xlabel('Multiples of Weight Vector')
ylabel('LRCs')
legend(sprintf('Dimension = %d',d))