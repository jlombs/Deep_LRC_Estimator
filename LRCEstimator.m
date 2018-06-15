function LRC = LRCEstimator(mu,nu,lambda,epsV,verb)
% Function which estimates the Littlewood-Richardson coefficient for
% weights mu,nu,lambda based on a randomized algorithm volume contraction
% algorithm. Weight vectors are double precision positive semidefinitie, and must be in weakly
% decreasing order while satisfying the summation rule. Verb is a flag for
% displaying progress data. eps is an accuracy parameter (between 0 and 1),
% with lower eps leading to longer but more accurate computations.

global findIndex bulkCoords n d eps

eps = epsV;

%% Establish variable shortcuts

n = numel(lambda); 
d = (n-1)*(n-2)/2; %Vector space dimension

buffer = ceil(10+d/10); % Autocorrelation buffer

%% Establish function for index finding for quick coordinate access given the data storage protocol

temp = reshape(1:(n+1)^2,[n+1,n+1])';
mask = flipud(tril(reshape(ones(1,(n+1)^2),[n+1,n+1])));
temp = temp.*mask;
indexes = find(temp);

findIndex = @(k,i) find(indexes == temp(k,i),1);

%% Generate matrix of bulk coordinates
% Each coordinate is an interior point in the LR-Hive, labeled by the
% column index k and row index i, with the origin at the bottom left
% corner.

numBulks = (n-1)*(n-2)/2;
bulkCoords = zeros(numBulks,2);
counter = 0;
for i = 1:(n-2)
   
    bulkCoords((counter+1):(counter + i),:) = [(1:i)',ones(i,1)*((n-1)-i)]+1;
    counter = counter + i;
    
end

bulkInds = zeros(numBulks,1);

for aa = 1:numBulks

    k = bulkCoords(aa,1);
    i = bulkCoords(aa,2);

    bulkInds(aa) = findIndex(k,i);
    
end

%% Preliminary check to see if the boundaries can form a hive

if verb
    
    fprintf('Performing Boundary Checks \n')
    
end
if sum(mu) + sum(nu) ~= sum(lambda)

    warning('Boundary data does not satisfy the Horn equality')
    return
end

%% Precheck on boundary format

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

%% Compute initial hive from the LP 

if verb
    
    fprintf('Testing for initial hive \n')
    
end

[~,failureFlag,LP] = generateMaxHive(mu,nu,lambda,0);

LPs = LP;

if failureFlag

    LRC = 0;
    warning('\n Maximal LP was not solvable--there may be no integer points within the hive boundaries')
    
    return
    
end

%% Progressively constrain hive inequalities until no solution is found, and store LP variables 

if verb
    
    fprintf('Finding hive contractions \n')
    
end

shift = 0;
while ~failureFlag
    
    shift = shift - 1;
    [~,failureFlag,LP] = generateMaxHive(mu,nu,lambda,shift);
    
    LPs = [LPs,LP]; %#ok

end

%% Enlarge to non-empty hive 
% (as a function of how many shifts--at least first nontrivial, but want larger to handle cut off corners etc...but larger is more time intensive as one should walk
% longer to ensure a better estimate of the number of integer points in this contracted hive, so a heuristic balance must be struck...)

shift = min([0,ceil(shift+2)]);

%% Plot testing
% Uncomment to look at the hive polytopes for bulk coordinates that live in
% R^3. (weight vectors having length 4)
%{
if size(LPs(1).A,2) == 3

    [V,~] = con2vert(LPs(1).A,LPs(1).b);
    k=convhulln(V);
    figure
    hold on
    for i=1:length(k)

        patch(V(k(i,:),1),V(k(i,:),2),V(k(i,:),3),'k','FaceAlpha',.5)

    end

    [V,~] = con2vert(LPs(abs(shift)+1).A,LPs(abs(shift)+1).b);
    kp=convhulln(V);
    hold on
    for i=1:length(kp)

        patch(V(kp(i,:),1),V(kp(i,:),2),V(kp(i,:),3),'r','FaceAlpha',.5)

    end

    title('Contracted and Original Interior Hive Polytope')
    
end
%}

if shift == 0
    
    %% Hive Small Enough that a Unique-Accumulating walk is Most Feasible
    
    if verb
        
        fprintf('Performing unique-accumulating computation due to hive size\n')
        
    end
    tmp = hiveWalkUniques(mu,nu,lambda,shift,inf,buffer,true,[]);
    LRC = size(tmp,2);
    
else
       
    %% Count the number of unique integer points one can find in the contracted volume using coordinate hit and run

    if verb
        
        fprintf('Approximating the lattice volume of the smallest hive \n')
        
    end

    [tmp,~] = hiveWalkUniques(mu,nu,lambda,shift,inf,buffer,true,[]);
    LRC = size(tmp,2);

    %% Set the random walk lengths for the next volume estimation (heuristic based on Ge,Ma and Zhang)

    N = 1600 * abs(shift);

    %% Now find all volume ratios for each enlarged polytop up to the original

    % Perform largest random walk in P0, the original polytope
    
    if verb
        
        fprintf('Walking through the original hive \n')
        
    end

    hives = hiveWalkUniques(mu,nu,lambda,0,N,buffer,false,tmp(:,randi(LRC)));

    % Since the polytopes are dilations, we compute every other
    % contraction ratio instead of each contraction by 1, reducing
    % accumulated errors while maintaining that the volume ratios are not
    % too extreme which would lead to inaccuracies in the sampling for
    % high dimensional polytopes
    
    x = -2;

    for i = -1:x:shift
        
        if verb
            
            fprintf('Estimating Volume Ratio %d/%d \n',abs(i),abs(shift))
            
        end
        
        % Find points inside new shell
        A = LPs(abs(i)+1).A * LPs(abs(i)+1).N'; 
        b = LPs(abs(i)+1).b;
                        
        ptInside = false(size(hives,2),1);
        
        for xx = 1:size(hives,2)
            
            ptInside(xx) = all(A *hives(:,xx) <= b);
            
        end
         
        % Update with volume ratio
        tot = sum(ptInside);
        LRC = LRC * size(hives,2)/tot;
 
        % Early exit if at the end
        
        if i == shift

            break

        end

        % Replace all outside points with new points randomly sampled inside the next shell
        % for the next iteration. This preservation of points and the process of starting at the largest
        % and working inward saves us computational time.
                
        if N - tot > 0
            
            outsideList = find(~ptInside);     
            tmp = find(ptInside,1);
            %tmp = tmp(randi(numel(tmp)));
            Y = hiveWalkUniques(mu,nu,lambda,i,N-tot,buffer,false,hives(:,tmp));
            hives(:,outsideList(1:N-tot)) = Y;
            hives(:,outsideList((N-tot+1):end)) = [];
            
        else
            
            hives(:,~ptInside) = [];
            
        end
        
        % Handle the last case if it would have been skipped through the x
        % jumps
        if i+x < shift
            
            if verb
                
                fprintf('Computing Volume Ratio %d/%d \n',abs(shift),abs(shift))
                
            end

            A = LPs(abs(shift)+1).A* LPs(abs(shift)+1).N';
            b = LPs(abs(shift)+1).b;
            ptInside = false(size(hives,2),1);
            
            for xx = 1:size(hives,2)
            
                ptInside(xx) = all(A *hives(:,xx) <= b);
            
            end
        
            LRC = LRC * size(hives,2)/sum(ptInside);

            break

        end

    end
    
end

LRC = round(LRC);

end

function [hives,exitStep] = hiveWalkUniques(mu,nu,lambda,shift,steps,buffer,uniqueFlag,hiveInit)
% Function which performs a coordinate hit and run walk on the hive
% polytope perscribed by the weight vectors mu nu and lambda (assumed to
% have the proper form already, see calling functions). Shift is a double
% specifying the polytope contraction, steps is a double specifying the
% number of c.h.a.r steps for a uniqueness run, or the number of samples
% desired in a blind run. buffer is a double specifying autocorrelation
% pause in number of steps between saved states for a blind run. uniqueFlag
% is a logical to trigger gathering only unique elements based on steps vs
% gather samples which may not be unique. hiveInit is an initial
% hive to start the walk at, empty if not provided. See calling functions for
% form of such an initial hive. 

global findIndex bulkCoords n eps 

exitStep = 1;

%% Establish variable shortcuts

totElem = (n+2)*(n+1)/2;
numBulks = size(bulkCoords,1);

%% Compute initial hive from the integer LP or use provided hive (one is guaranteed to exist by the calling structure)

initLength = 100;
hives = zeros(totElem,initLength);
hiveNum = 1;
    
if isempty(hiveInit)

    [Hijk,~] = generateMaxHive(mu,nu,lambda,shift);
    
else
    
    Hijk = hiveInit;
    
end

hives(:,1) = Hijk;

%% Check for a tight hive and also store bulk indices and the neighbor indices of each bulk for later reference
% A tight hive means that there are no flexible coordinates in the polytope
% bulk, which either indicates the existance of only 1 integer point, or
% there could be others which are not axis aligned. This is discussed in
% the readme paper as a shortcoming of this bulk-index c.h.a.r. implementation. 
   
bulkInds = zeros(numBulks,1);
neighborStore = zeros(numBulks,12);

data.tight = true;

for aa = 1:numBulks

    k = bulkCoords(aa,1);
    i = bulkCoords(aa,2);

    bulkInds(aa) = findIndex(k,i);
    
    neighborStore(aa,1:9)  = [findIndex(k-1,i+1),findIndex(k,i+1),findIndex(k-1,i+2),findIndex(k+1,i),findIndex(k+1,i-1),...
        findIndex(k+2,i-1),findIndex(k-1,i),findIndex(k,i-1),findIndex(k-1,i-1)];
    
    %See commented section at the end of this function for computational
    %details: this is the quicker version but more obtuse
    
    a = Hijk(neighborStore(aa,1));
    b = Hijk(neighborStore(aa,2));
    c = Hijk(neighborStore(aa,3));
    d = Hijk(neighborStore(aa,4));
    e = Hijk(neighborStore(aa,5));
    f = Hijk(neighborStore(aa,6));
    g = Hijk(neighborStore(aa,7));
    h = Hijk(neighborStore(aa,8));
    j = Hijk(neighborStore(aa,9));
    
    maxBound =  min([a + b - c,...
        d + e - f,g + h - j,...
        ]) + shift;
    minBound = max([a + h - g,...
        e + b - d, a + d - b,...
        h + d - e,g + b - a,...
        g + e - h]) - shift;
    
    if i>2
        
        neighborStore(aa,10) = findIndex(k+1,i-2);
        
    end
    
    if k+i~=n+1
        
        neighborStore(aa,11) = findIndex(k+1,i+1);
        
    end
    
    if k>2
       
        neighborStore(aa,end) = findIndex(k-2,i+1);
        
    end

    if maxBound ~= minBound
        
        if i > 2

                maxBound = min([maxBound, h + e - Hijk(neighborStore(aa,10)) + shift]);

        end

        if maxBound ~= minBound

            if k+i ~= n+1
           
                maxBound = min([maxBound, b + d - Hijk(neighborStore(aa,11)) + shift]);
            
            end

            if maxBound ~= minBound

                if k > 2

                    maxBound = min([maxBound, a + g -  Hijk(neighborStore(aa,end)) + shift]);           

                end

            end
            
        end
        
    end
    
    % No tight hive, establish starting coordinate and set flag
    if minBound ~= maxBound

        coord = aa;
        data.tight = false;

    end

    % Tight hive, exit
    if aa == numBulks && data.tight == true
        
        hives = hives(:,1);
        return

    end

end

%% Start the walk

stepCount = 0;
bufferCount = 1;

if uniqueFlag
    
    uniqueBootO = 1;
    
end

while hiveNum < steps
    
    %% Pick first flexible instance as found from tight check as initial condition, otherwise randomly choose a new coordinate

    stepCount = stepCount + 1;
    
    if stepCount > 1
        
        coord = randi(numBulks);
        coordIndex = bulkCoords(coord,:);
        k = coordIndex(1);
        i = coordIndex(2);
    
        %% Find limits of polytope along this direction

        a = Hijk(neighborStore(coord,1));
        b = Hijk(neighborStore(coord,2));
        c = Hijk(neighborStore(coord,3));
        d = Hijk(neighborStore(coord,4));
        e = Hijk(neighborStore(coord,5));
        f = Hijk(neighborStore(coord,6));
        g = Hijk(neighborStore(coord,7));
        h = Hijk(neighborStore(coord,8));
        j = Hijk(neighborStore(coord,9));
        
        
        maxBound =  min([a + b - c,...
            d + e - f,g + h - j,...
            ]) + shift;
        minBound = max([a + h - g,...
            e + b - d, a + d - b,...
            h + d - e,g + b - a,...
            g + e - h]) - shift;

        if minBound ~= maxBound

            if i > 2

                maxBound = min([maxBound, h + e - Hijk(neighborStore(coord,10)) + shift]);

            end

            if maxBound ~= minBound

                if k+i ~= n+1

                    maxBound = min([maxBound, b + d - Hijk(neighborStore(coord,11)) + shift]);
                    
                end

                if maxBound ~= minBound

                    if k > 2
                        
                        maxBound = min([maxBound, a + g - Hijk(neighborStore(coord,12)) + shift]);           

                    end

                else
                    
                    continue
                    
                end
                
            else
                
                continue
                
            end
            
        else
            
            continue
            
        end

        if minBound == maxBound

            continue

        end

    end
    
    %% Perturb uniformly randomly within the limits
    
    Hijk(bulkInds(coord)) = randi([minBound,maxBound]);
    
    %% Check if the hive is unique
    
    if uniqueFlag
        
            %% Potential early break if not finding new points after ceil(steps/buffer) number of attempts
            % Note, full rank hive polytopes are not excessively narrow,
            % and as a result, this heuristic seems ok. If one had to worry
            % about wandering through a bottleneck, the full sample time
            % would be necessary.          

        if stepCount - uniqueBootO > ceil(buffer/eps^2)

            exitStep = stepCount - ceil(buffer/eps^2);
            break

        end
        
        if ~ismemberR2012aR(Hijk,hives(:,1:hiveNum))
            
            uniqueBootO = stepCount;

            %% Store unique data
            
            hiveNum = hiveNum + 1;

            if hiveNum > size(hives,2)

                hives = padarray(hives,[0,initLength],0,'post');  
            end

            hives(:,hiveNum) = Hijk;    
                        
        end
           
    else 
        
        bufferCount = bufferCount + 1;
        
        if bufferCount == buffer
            
            bufferCount = 1;
            hiveNum = hiveNum + 1;

            if hiveNum > size(hives,2)

                hives = padarray(hives,[0,initLength],0,'post');  
            end

            hives(:,hiveNum) = Hijk;
            
        end
        
    end

end

hives = hives(:,1:hiveNum);

%{

%% Detailed Discussion of what's happening above: combining all maxmin
tests into vectorized form with quick stored variables

% Four vertical rhombuses

    %Up

        %All bulk coeffs admit an up rhombus
        %Hijk(k,i) + Hijk(k-1,i+2) <= Hijk(k-1,i+1) + Hijk(k,i+1)
        maxBound = Hijk(k-1,i+1) + Hijk(k,i+1) - Hijk(k-1,i+2)+shift;

    %Left

        %All bulk coeffs admit a left rhombus
        %Hijk(k-1,i+1) + Hijk(k,i-1) <= Hijk(k,i) + Hijk(k-1,i)
        minBound = Hijk(k-1,i+1) + Hijk(k,i-1) - Hijk(k-1,i)-shift;

    %Early boot

    if maxBound == minBound

        return

    end

    %Down

        %If the bottom exists, then all exist and can check
        if i > 2

            %Hijk(k,i) + Hijk(k+1,i-2) <= Hijk(k,i-1) + Hijk(k+1,i-1)
            maxBound = min([maxBound, Hijk(k,i-1) + Hijk(k+1,i-1) - Hijk(k+1,i-2)+shift]);
            
            %Early boot

            if maxBound == minBound

                return

            end

        end

    %Right

        %All bulk coeffs admit a right rhombus
        %Hijk(k+1,i-1) + Hijk(k,i+1) <= Hijk(k,i) + Hijk(k+1,i)
        minBound = max([minBound, Hijk(k+1,i-1) + Hijk(k,i+1) - Hijk(k+1,i)-shift]);

    %Early boot

    if maxBound == minBound

        return

    end

% Four acute rhombuses

    %Right Up

        %All bulk except for the right interior column have right-up
        if k+i ~= n+1

            %Hijk(k,i) + Hijk(k+1,i+1) <= Hijk(k,i+1) + Hijk(k+1,i)
            maxBound = min([maxBound, Hijk(k,i+1) + Hijk(k+1,i) - Hijk(k+1,i+1)+shift]);
        
            %Early boot

            if maxBound == minBound

                return

            end
            
        end

    %Right Down

        %All bulk admit a right-down rhombus
        %Hijk(k,i) + Hijk(k+2,i-1) <= Hijk(k+1,i) + Hijk(k+1,i-1)
        maxBound = min([maxBound, Hijk(k+1,i) + Hijk(k+1,i-1) - Hijk(k+2,i-1)+shift]);
        
    %Early boot

    if maxBound == minBound

        return

    end

    %Left Up

        %All bulk except for the left interior column have left-up
        if k > 2

            %Hijk(k,i) + Hijk(k-2,i+1) <= Hijk(k-1,i+1) + Hijk(k-1,i)
            maxBound = min([maxBound, Hijk(k-1,i+1) + Hijk(k-1,i) - Hijk(k-2,i+1)+shift]);           
        
            %Early boot

            if maxBound == minBound

                return

            end
            
        end

    %Left Down

        %All bulk admit a left-down
        %Hijk(k,i) + Hijk(k-1,i-1) <= Hijk(k-1,i) + Hijk(k,i-1)
        maxBound = min([maxBound, Hijk(k-1,i) + Hijk(k,i-1) - Hijk(k-1,i-1)+shift]);     

    %Early boot

    if maxBound == minBound

        return

    end

% Four obtuse rhombuses

    %Right Up

        %All bulk admit a right-up
        %Hijk(k-1,i+1) + Hijk(k+1,i) <= Hijk(k,i) + Hijk(k,i+1)
        minBound = max([minBound, Hijk(k-1,i+1) + Hijk(k+1,i) - Hijk(k,i+1)-shift]);
        
    %Early boot

    if maxBound == minBound

        return

    end 

    %Right Down

        %All bulk admit a right-down
        %Hijk(k,i-1) + Hijk(k+1,i) <= Hijk(k,i) + Hijk(k+1,i-1)
        minBound = max([minBound, Hijk(k,i-1) + Hijk(k+1,i) - Hijk(k+1,i-1)-shift]);
        
    %Early boot

    if maxBound == minBound

        return

    end

    %Left Up

        %All bulk admit a left-up
        %Hijk(k-1,i) + Hijk(k,i+1) <= Hijk(k,i) + Hijk(k-1,i+1)
        minBound = max([minBound, Hijk(k-1,i) + Hijk(k,i+1) - Hijk(k-1,i+1)-shift]);
        
    %Early boot

    if maxBound == minBound

        return

    end

    %Left Down

        %All bulk admit a left-down
        %Hijk(k-1,i) + Hijk(k+1,i-1) <= Hijk(k,i) + Hijk(k,i-1)
        minBound = max([minBound, Hijk(k-1,i) + Hijk(k+1,i-1) - Hijk(k,i-1)-shift]);

 
end
%}
end

function [lia] = ismemberR2012aR(a,b)
% Shortcut membership test with less overhead, given our known inputs

uA = repmat(a,1,size(b,2));
lia = any(all(uA==b,1));

end


function [A,b,AEQ,bEQ] = hiveLP(mu,nu,lambda,s)
% Generates the Hive LP given input boundary vectors.
% Assumes vectors are pre-sorted integers in weakly decreasing order.
% s is a shift on the polytope inequalities, s=0 for original hive.

global n findIndex

%% Shortcut variables
totElem = (n+1)*(n+2)/2;

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

function [x,failureFlag,P] = generateMaxHive(mu,nu,lambda,shift)
%Function which solves the max (integer by theory) LP for boundary vectors
%mu,nu,lambda assumed to be already in weakly decreasing interger order.
%Shift is a double parameter which gives the polytope contraction. Outputs
%are x as the solution, failureFlag if no solution was found, and P as an
%LP structure with fields A and b.

global n

failureFlag = 0;

%% Build out the LP structure

[P.A,P.b,A_eq,b_eq] = hiveLP(mu,nu,lambda,shift);


%% Perform the linear program and round the answer (as we know analytically the results should be integers)
options = optimoptions(@linprog,'Display','off','Algorithm','dual-simplex');
x = linprog(-ones((n+1)*(n+2)/2,1),P.A,P.b,A_eq,b_eq,[],[],options);

if isempty(x)
    
    failureFlag = 1;
    x = [];
    P = [];
    return
    
else

    x = round(x);
    
end

%% Project LP into minimal form and save for output
N = null(A_eq);
z = linsolve(A_eq, b_eq);
P.b = P.b - P.A * z;
P.A = P.A*N;
P.N = N;

end