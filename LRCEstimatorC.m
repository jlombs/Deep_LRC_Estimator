function LRC = LRCEstimatorC(mu,nu,lambda,eps,verb)
% Function which computes the LRC based on 1-step volume estimations using the
% Cousins algorithm for weight vectors mu, nu, and lambda with precision
% eps [0,1] (lower eps, higher precision). verb is a verbosity flag for
% console output.


%% Establish variable shortcuts

n = length(lambda); 

d = (n-1)*(n-2)/2; %Vector space dimension

if d == 1
    
    x = log(2);
    
else
    
    x = log(d);
    
end

%% Preliminary check to see if the boundaries can form a hive

if verb
    
    fprintf('Performing Boundary Checks \n')
    
end

if numel(nu) > n || numel(mu) > n

    error('lambda must be the largest vector, or all vectors must be the same length')

end

if size(nu,2) < size(nu,1)
    
    nu = nu';
    
end

if size(mu,2) < size(mu,1)
    
    mu = mu';
    
end

if size(lambda,2) < size(lambda,1)
    
    lambda = lambda';
    
end

if numel(nu) < n

    nu = padarray(nu,[0,n-numel(nu)],0,'post');

end

if numel(mu) < n

    mu = padarray(mu,[0,n-numel(mu)],0,'post');

end

if ~isa(mu,'double') || ~isa(nu,'double') || ~isa(lambda,'double')

    error('mu, nu, and lambda must all be integer arrays of double datatype')

end

if any(mod(mu,1)) || any(mod(nu,1)) || any(mod(lambda,1))  

    error('mu, nu, and lambda must all be integer arrays of double datatype')

end

mu = sort(mu,'descend');
nu = sort(nu,'descend');
lambda = sort(lambda,'descend');

if mu(end) < 0 || nu(end) < 0 || lambda(end) < 0 || mu(1) == 0 || nu(1) == 0 || lambda(1) == 0

    error('mu, nu, and lambda must be positive vectors with at least 1 nonzero term')

end   

if sum(mu) + sum(nu) ~= sum(lambda)
    
    error('mu, nu, and lambda must satisfy the Horn equality')
    
end

%% Compute Minimal LPs for the enlarged polytope Q and the original hive P

if verb
    
    fprintf('Computing polytopes and processing \n')
    
end

[Q.A,Q.b,Q.A_eq,Q.b_eq] = hiveLP(mu,nu,lambda,1);
[P.A,P.b,P.A_eq,P.b_eq] = hiveLP(mu,nu,lambda,0);

%% Only project out the equality subspace, and don't look for 0-width facets as there shouldn't be any
% Note, this code is partially reproduced from Cousins 'preprocess.m'
% script. We treat the interior point finding differently and avoid
% centering the polytopes by translation (as is done in the preprocess script by default) as this is not a lattice volume
% preserving transformation.

N = null(P.A_eq);
z = linsolve(P.A_eq, P.b_eq);
Pp.b = P.b - P.A * z;
Pp.A = P.A*N;

N = null(Q.A_eq);
z = linsolve(Q.A_eq, Q.b_eq);
Qp.b = Q.b - Q.A * z;
Qp.A = Q.A*N;

Q = Qp;
P = Pp;

if verb

    fprintf('Trying to find a point inside the convex body...\n');

end

dim = size(Q.A,2);
p=zeros(dim,1);
opts = optimoptions('linprog','Algorithm','dual-simplex','Display','off');

for i=1:size(Q.A,2)

    [~,f_index]=min(Q.b-Q.A*p);
    f = Q.A(f_index,:);
    [y]=linprog(f,Q.A, Q.b,[],[],[],[],[],opts);
    p = ((i-1)*p+y)/i;

    if verb

        if mod(i,10)==0
            fprintf('%d its, %f frac of equations satisfied (this one has %f)\n', i, sum(Q.A*p<=(Q.b))/length(Q.b),sum(Q.A*y<=(Q.b))/length(Q.b));
        end

    end

end

Q.p = p;

%% Plot testing
% Uncomment to look at the hive polytopes for bulk coordinates that live in
% R^3. (weight vectors having length 4)

%{
if size(P.A,2) == 3

    [V,~] = con2vert(Q.A,Q.b);
    k=convhulln(V);
    figure
    hold on
    for i=1:length(k)

        patch(V(k(i,:),1),V(k(i,:),2),V(k(i,:),3),'k','FaceAlpha',.5)

    end

    [V,~] = con2vert(P.A,P.b);
    kp=convhulln(V);
    hold on
    for i=1:length(kp)

        patch(V(kp(i,:),1),V(kp(i,:),2),V(kp(i,:),3),'r','FaceAlpha',.5)

    end

    title('Expanded and Original Interior Hive Polytope')
    
end

%}

%% Compute the larger volume (Q) with a coordinate random walk

if verb 
    
    fprintf('Computing the polytope volume \n')
    
end

volFlags = sprintf('-round -walk char -verb %d',verb);
hiveVolume = Volume(Q,[],eps,volFlags);

%% Compute fraction of lattice points inside P from a sampling on Q

if verb 
    
    fprintf('Computing the point ratio \n')
    
end

%Set variables for fractional convergence within eps and initialize vars

err = inf;
fracO = inf;

sampleChunk = 1000;
numInside = 0;
totSamples = 0;

%Adaptive centering hit and run, invariant under linear transforms
options.method = 'achr';

while err > eps
   
    %Take sampleChunk samples at a time from Q, picking up where you leave off
    %(unless its the first sample, which starts at the chebycenter and goes
    %from there in the uniform sampler)
    
    if err == inf
        
        X = round(cprnd(sampleChunk,Q.A,Q.b))';
        
    else
        
        options.x0 = X(:,end);
        X = round(cprnd(sampleChunk,Q.A,Q.b,options))';
        
    end
    
    totSamples = totSamples + sampleChunk;
    
    %Compute the number of points inside P
    numInside = numInside + sum(arrayfun(@(x) all(P.A*X(:,x)<=P.b),1:sampleChunk)); 
    
    %Update the fractions and error measure
    frac = numInside/totSamples;
    
    err = abs(fracO - frac)/frac;
    
    fracO = frac;
    
end

%% Output the LRC

LRC = round(hiveVolume*frac);

end

function [A,b,AEQ,bEQ] = hiveLP(mu,nu,lambda,s)
%Generates the LP given input boundary vectors.
% Assumes vectors are pre-sorted integers in weakly decreasing order.
% s is a shift on the polytope inequalities, s=0 for original hive.

%% Shortcut variables

n = length(mu);
totElem = (n+1)*(n+2)/2;

%% Establish function for index finding for quick coordinate access given the data storage protocol

temp = reshape(1:(n+1)^2,[n+1,n+1])';
mask = flipud(tril(reshape(ones(1,(n+1)^2),[n+1,n+1])));
temp = temp.*mask;
indexes = find(temp);

findIndex = @(k,i) find(indexes == temp(k,i),1);

%% Create storage

A = zeros(n*(n-1)*3/2,totElem);
b = zeros(n*(n-1)*3/2,1)+s;

AEQ = zeros(3*n,totElem);
bEQ = zeros(3*n,1);

%% Build boundary values

csMu = cumsum(mu);
csNu = cumsum(nu);
csLambda = cumsum(lambda);

%% Loop over boundary and interior to build out the LP

counterB = 1;
counterI = 1;

for k = 0:n
    
    for i = 0:(n-k)
        
        % Boundary
        if k == 0 && i == 0
            
            AEQ(counterB,1) = 1;
            counterB = counterB + 1;
            
        elseif k == 0 && i ~= 0
            
            AEQ(counterB,findIndex(k+1,i+1)) = 1;
            bEQ(counterB) = csMu(i);
            counterB = counterB + 1;
            
        elseif i == 0 && k ~= 0
            
            AEQ(counterB,findIndex(k+1,i+1)) = 1;
            bEQ(counterB) = csLambda(k);
            counterB = counterB + 1;
            
        elseif k == n-i
            
            AEQ(counterB,findIndex(k+1,i+1)) = 1;
            bEQ(counterB) = csMu(end) + csNu(k);
            counterB = counterB + 1;
            
        end
        
        % Interior 
        
        %Right Hijk(k+1,i+1) ~= 0 && Hijk(k,i+1) ~= 0 && Hijk(k+1,i) ~=0 && Hjk + Hijk(k+1,i+1) > Hijk(k,i+1)+Hijk(k+1,i)
        if i+k < n-1

            A(counterI,findIndex(k+1,i+1)) = 1;
            A(counterI,findIndex(k+2,i+2)) = 1;
            A(counterI,findIndex(k+1,i+2)) = -1;
            A(counterI,findIndex(k+2,i+1)) = -1;
            
            counterI = counterI + 1;
            

        end
        %Top Hijk(k+1,i-2) ~= 0 && Hijk(k,i-1) ~= 0 && Hijk(k+1,i-1) ~=0 && Hjk + Hijk(k+1,i-2) > Hijk(k,i-1)+Hijk(k+1,i-1)
        if i > 1

            A(counterI,findIndex(k+1,i+1)) = 1;
            A(counterI,findIndex(k+2,i-1)) = 1;
            A(counterI,findIndex(k+1,i-0)) = -1;
            A(counterI,findIndex(k+2,i-0)) = -1;
            
            counterI = counterI + 1;

        end
        %Left Hijk(k-2,i+1) ~= 0 && Hijk(k-1,i) ~= 0 && Hijk(k-1,i+1) ~=0 && Hjk + Hijk(k-2,i+1) > Hijk(k-1,i)+Hijk(k-1,i+1)
        if k > 1

            A(counterI,findIndex(k+1,i+1)) = 1;
            A(counterI,findIndex(k-1,i+2)) = 1;
            A(counterI,findIndex(k-0,i+1)) = -1;
            A(counterI,findIndex(k-0,i+2)) = -1;
            
            counterI = counterI + 1;

        end
        
    end

end

end