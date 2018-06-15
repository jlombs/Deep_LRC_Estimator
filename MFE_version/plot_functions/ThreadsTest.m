function [volumes] = ThreadsTest(P,E,eps,flags,num_steps, step_size)

%---INPUT VALUES---
%P: the polytope, {x | P.A*x <= P.b}
%E: the ellipsoid, {x | (x-E.v)'E.E(x-E.v)<=1}
%eps: the target relative error
%flags: a string of input flags. see parseFlags.m
%num_steps: the number of steps to allocate per phase
%step_size: the size of each step

%---RETURN VALUES---
%volumes: a num_steps X 1 array where volumes(i) will be the volume after
%         running i*step_size steps per phase.

%assign default values if not assigned in function call
if exist('flags','var')==0
    flags = '';
end
if exist('eps','var')==0
    eps = 0.20;
end
if exist('E','var')==0
    E=[];
end

%prepare our objects
K = ConvexBody(P,E,eps,flags);
if K.verb>=1
    fprintf('--------%d-Dimension Convex Body------\n\n', K.dim);
end

%make sure the provided point is inside the body
if ~in_K(K,zeros(K.dim,1))
    error('The point provided is not in the convex body! Please provide a different center point.');
end

%let's initialiez some helpful constants
[ratio,num_threads,C] = assignConstants(K);

if isKey(K.flagmap,'round')
    %rounding phase
    
    %round the body once as a preprocessing step
    %note that K is modified inside round()
    [T]=round(K,num_threads);
else
    %we are not rounding the body
    T=eye(K.dim);
end

%make it so the min distance to dK from 0 is exactly 1
%helpful for sampler
% [det_round]=unitBallify(K,det_round);

if K.verb>=1
    fprintf('------Volume Start------\n');
end

%these are our starting points for each thread
%we shifted our body so that K contains the origin, so x \in K
x = zeros(K.dim,num_threads);

%compute the annealing schedule that keeps E(Y^2)/E(Y)^2<=C
[a_sched] = getAnnealingSchedule(K,ratio,num_threads,C);
K.m = length(a_sched);

%compute the initial volume, multiplied by the determinant of the rounding
% matrix.
volume = (pi/a_sched(1))^(K.dim/2)*abs(det(T));

%initialize additional helpful variables
volumes=volume*ones(num_steps,1);

if K.verb>=1
    fprintf('Num Phases: %d\n', K.m);
end

for i=1:length(a_sched)-1
    if K.verb>=1
        fprintf('Phase %d Volume: %e', i-1,volumes(end));
        if K.verb>=2
            %compute how many "ratios" we've stepped down so far
            sigma_index=round(log(a_sched(i)/a_sched(1))/log(ratio));
            fprintf(',   sigma_%d=%e', sigma_index,1/sqrt(2*a_sched(i)));
        end
        fprintf('\n');
    end
    
    %take the number of steps in this phase
    
    [ratios,x]=MixStepsPhase(K,x,a_sched(i+1),a_sched(i),num_steps, step_size);
    volumes=volumes.*ratios;
    
end

if K.verb>=1
    fprintf('------Volume Complete------\n\n');
    fprintf('Final Volume: %e,  final sigma=%e\n',volumes(end), 1/sqrt(2*a_sched(end)));

end

end

function [ratio,num_threads,C] = assignConstants(K)
%initialize some hard-coded constants
if isKey(K.flagmap,'ratio')
    ratio = K.flagmap('ratio');
else
    ratio = 1-1/K.dim;
end
if isKey(K.flagmap,'num_t')
    num_threads = K.flagmap('num_t');
else
    num_threads = 5;
end
if isKey(K.flagmap,'C')
    C = K.flagmap('C');
else
    C = 2;
end
end