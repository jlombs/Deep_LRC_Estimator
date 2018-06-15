function [ P ] = preprocess( P ,verbFlag)%, to_round)
%PREPROCESS Takes in an object P, which is a polytope with the following
%properties:
%
%the polytope is {x | P.A * x <= P.b} \cap {x | P.A_eq * x = P.b_eq}
%P.A
%P.b
%P.A_eq
%P.b_eq
%
%an optional argument P.p, which is a point satisfying P.A * P.p <= P.b
%
% The algorithm will first restrict to the equality subspace defined by
% P.A_eq, P.b_eq. Then, it will look at the width of the polytope in each
% facet direction. If it finds a width 0 facet, say with normal vector a_i,
% then the polytope actually lies in a lower dimensional subspace (defined
% by a_i and the point found by the LP solver when optimizing w.r.t. a_i).
% So, we further restrict based on all width 0 facets.
%
% After this preprocessing, you still may want to round the polytope for
% the volume algorithm to be accurate and efficient.
%
% Example usage:
%
% P = makeBody('cube',10);
% P = preprocess(P);
% vol = Volume(P);

%check to make sure P.A and P.b are defined, and appropriately sized

if (isfield(P,'A')==0 || isfield(P,'b')==0) || (isempty(P.A) || isempty(P.b))
    %either P.A or P.b do not exist
    error('You need to define both P.A and P.b for a polytope {x | P.A*x <= P.b}.');
end

num_constraints = size(P.A,1);
dim = size(P.A,2);

if verbFlag
    
    fprintf('Currently (P.A, P.b) are in %d dimensions\n', dim);
    
end

if size(P.b,2)~= 1 || num_constraints ~= size(P.b,1)
    error('Dimensions of P.b do not align with P.A.\nP.b should be a %d x 1 vector.',num_constraints);
end

if (isfield(P,'A_eq')==0 || isempty(P.A_eq)) && ...
        (isfield(P,'b_eq')==0 || isempty(P.b_eq))
    P.A_eq = [];
    P.b_eq = [];
else
    %we are in a smaller dimensional space.
    %let's do some linear algebra to restrict
    %Ax<=b to this space
    if verbFlag
        
        fprintf('Restricting {P.A, P.b} to the equality subspace.\n');
    
    end
    N = null(P.A_eq);
    z = linsolve(P.A_eq, P.b_eq);
    P.b = P.b - P.A * z;
    P.A = P.A*N;
    dim = size(P.A,2);
    
    if verbFlag
        
        fprintf('Now in %d dimensions after restricting.\n', dim);

    end
    
end

%check for width 0 facets to make sure we are full dimensional
%also check for feasibility

if verbFlag
    
    fprintf('Checking for width 0 facets...\n');

end

[widths, vals] = getWidths(P,verbFlag);
eps_cutoff = 1e-7;

num_eq = sum(widths<eps_cutoff);

if num_eq > 0
    
    if verbFlag
        
        fprintf('Found %d width 0 facets.\n', num_eq);
    
    end
    
    eq_constraints = zeros(num_eq,dim);
    eq_rhs = zeros(num_eq,1);
    curr = 1;
    for i=1:length(widths)
        if abs(widths(i)) < eps_cutoff
            %this facet has width 0
            %we add this constraint to the equality
            %subspace, with the value x_opts(i) defining
            %the shift of the facet
            
            eq_constraints(curr,:) = P.A(i,:);
            eq_rhs(curr) = vals(i);
            curr = curr + 1;
        end
    end
    %restrict to the degenerate subspace
    N = null(eq_constraints);
    z = linsolve(eq_constraints, eq_rhs);
    P.b = P.b - P.A * z;
    P.A = P.A*N;
end

%remove zero rows
new_A = [];
new_b = [];
for i=1:size(P.A,1)
   if norm(P.A(i,:))>1e-10
       new_A = [new_A; P.A(i,:)];
       new_b = [new_b; P.b(i)];
   end
end

if verbFlag
    
    fprintf('Removed %d zero rows\n', size(P.A,1)-size(new_A,1));
    
end

P.A = new_A;
P.b = new_b;
dim = size(P.A,2);

if verbFlag

    fprintf('Trying to find a point inside the convex body...\n');
    
end

% if isfield(P,'p')==0 || isempty(P.p)
p=zeros(dim,1);
%options = optimset('Display','none');
options = optimoptions('linprog','Algorithm','dual-simplex','Display','off');

for i=1:size(P.A,2)
%     f = randn(dim,1);
    [~,f_index]=min(P.b-P.A*p);
    f = P.A(f_index,:);
    [y]=linprog(f,P.A, P.b,[],[],[],[],[],options);
    p = ((i-1)*p+y)/i;
    
    if verbFlag
        
        if mod(i,10)==0
            fprintf('%d its, %f frac of equations satisfied (this one has %f)\n', i, sum(P.A*p<=(P.b))/length(P.b),sum(P.A*y<=(P.b))/length(P.b));
        end
        
    end
    %     fprintf('Number of low rows: %d\n', sum(P.b - P.A*p < eps_cutoff));
end

%hopefully found a point reasonably inside the polytope
%shift so that this point is the origin
P.b = P.b - P.A*p;

if min(P.b) < -eps_cutoff
   error('We tried to find a point inside the polytope but failed.'); 
elseif min(P.b) < 0
    
    if verbFlag
        
        fprintf('The point is very close to inside, so we slightly perturb P.b to make it be so.\n');
        
    end
    P.b = P.b + eps_cutoff;
else
    
    if verbFlag
        
        fprintf('We found a point inside the polytope!\n');
        
    end
end

if verbFlag
    
    fprintf('Final polytope is in %d dimensions.\n', dim);
    
end

%ORIGINAL CODE WAS MISSING THIS LINE
%P.p = p;

end

function [widths, vals] = getWidths(P,verbFlag)

%options = optimset('Display','none');
options = optimoptions('linprog','Algorithm','dual-simplex','Display','off');

widths = zeros(size(P.A,1),1);
vals = zeros(size(P.A,1),1);

if verbFlag
    
    fprintf('%d constraints to check.\n', length(widths));
    
end

chunk = ceil(length(widths)/10);

for i=1:length(widths)
    %[~,max_dist] = linprog(f, [A1; K.A(i,:) -1],[K.b; K.b(i)],[],[],[],[],[],options);
    %[~,min_dist] = linprog(f, [A1; K.A(i,:) 1],[K.b; K.b(i)],[],[],[],[],[],options);
    if verbFlag
        
        if mod(i,chunk)==0
            fprintf('%d/%d done..', i, length(widths));
        end
        
    end

    [x,max_dist,exitflag] = linprog(-P.A(i,:), P.A, P.b,[],[],[],[],[],options);
    
    if exitflag==-2
        error('The polytope provided has no feasible points! (so it appears)');
    end
    [y,min_dist] = linprog(P.A(i,:), P.A, P.b,[],[],[],[],[],options);
    
    
    widths(i) = abs(max_dist+min_dist)/norm(P.A(i,:),2);
    vals(i) = P.A(i,:)*x;
    %     if widths(i) < 1e-7
    %         fprintf('(%e,%e)\n', vals(i), min_dist);
    %     end
end

end