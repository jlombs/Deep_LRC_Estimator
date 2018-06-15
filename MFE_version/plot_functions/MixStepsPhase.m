function [ ratios, x ] = MixStepsPhase(K,x,a, last_a, num_steps, step_size)
%MixStepsPhase This function will take in a convex body K, current points
%for each thread x, next distribution a, distribution to sample from, and
%the number of steps to run for this phase

%it will return the ratio estimates for z(a)/z(last_a) at each step i for
%i=1, ..., num_steps.

[~,num_threads,~]=assignConstants(K);

ratios=zeros(num_steps,1);

fn = 0;

resetSlacks(K,x);

it=1;
for i=1:ceil(num_steps*step_size/num_threads)
    for j=1:num_threads
        if it>num_steps*step_size
            break;
        end
        x(:,j) = getNextPoint(K,x(:,j),last_a, j);
        
        fn = fn + eval_exp(x(:,j),a)/eval_exp(x(:,j),last_a);
        
        if mod(it,step_size)==0        
            ratios(it/step_size) = fn/it;
        end
        it = it+1;
        
    end
    
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

