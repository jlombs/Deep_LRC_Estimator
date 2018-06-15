%% Hive Integer Point Estimation Test Script

% This script is designed to test the computation of the
% Littlewood-Richardson coefficient for provided weight vectors mu, nu, and
% lambda. These row vectors must be positive semidefinite integers of double type in weakly
% descending order, with sum(mu)+sum(nu) = sum(lambda). The function called
% is a randomized algorithm based on a coordinate hit and run on a
% continuum version of the hive polytope. We use Ben Cousins algorithm for estimating the volume of an 
% enlarged hive, and compute the fractional number of lattice points to the original polytope using a rounded 
% adaptive centering hit and run.

clear
clc
close all

%% Boundary data

%{
mu = [40,30,20,10];
nu = mu;
lambda = [65,55,45,35];
% Note, should be 506

%}

%{
mu = [40,30,20,10];
nu = mu;
lambda = [65,55,46,34];
% Note, should be 505

%}

%{
mu = 7*[40,30,20,10];
nu = mu;
lambda = 7*[65,55,46,34];
% Note, should be 133792

%}

%{
mu = [2,1,1];
nu = mu;
lambda = [3,3,2];
n = 3;
% Note should be 1
%}

%{
mu = 10*[2,1,0];
nu = mu;
lambda = 10*[3,2,1];
n = 3;
% Note should be 11
%}

%{
mu = [42    39    31    29    17     1];
nu = [43    41    40    23     9     5 ];
lambda = [83    69    44    44    40    40];
% Note should be 156
%}

%{
mu = sort(randi(50,[1,6]),'descend');
nu = sort(randi(50,[1,6]),'descend');
lambda = sort(mu + nu(randperm(numel(nu))),'descend');
%}

%% Set tollerances and flags

eps = .5; % Eps is a fractional error parameter between 0 and 1
verbosityFlag = true; % verbosityFlag indicates console feedback

LRC = LRCEstimatorC(mu,nu,lambda,eps,verbosityFlag);

fprintf('LRC is %d \n',LRC)
