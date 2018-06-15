%Script to test Linearized Least Squares fitting as an LP
close all

%% Generate random data to be fit
pts = 10;

x = sort(rand(1,pts));
f = @(x) 3*x.^2+9*x+6;
y = f(x);

%% Build LP for polyfit of order n
n = 2;
A = ones(pts,n+1);

for i = 1:pts
    A(i,2:end) = arrayfun(@(z) x(i)^z,1:n);
end
b = y;


%% Solve LP
[LP,resnorm,residual] = lsqnonneg(A,b');
LP = LP';

%% Plot
figure(1)
plot(x,y,'k-')
hold on

numSamp = 50;
xp = linspace(x(1),x(end),numSamp);
fp = ones(1,numSamp);
for i = 1:numSamp
    fp(i) = sum(LP.*arrayfun(@(z) xp(i)^z,0:n));
end
    
plot(xp,fp,'r--')

