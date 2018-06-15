%% Hive Integer Point Estimation Test Script

% This script is designed to test the computation of the
% Littlewood-Richardson coefficient for provided weight vectors mu, nu, and
% lambda. These row vectors must be positive semidefinite integers of double type in weakly
% descending order, with sum(mu)+sum(nu) = sum(lambda). The function called
% is a randomized algorithm based on a coordinate hit and run on the hive
% lattice itself, using telescoping lattice volume ratios and an inner
% computation of a small lattice hive found from polytope contractions.

clear
clc
close all

%% Boundary data

%%{
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

%% Set tollerances and flags

eps = .25; % Eps is a fractional error parameter between 0 and 1
verbosityFlag = true; % verbosityFlag indicates console feedback

%% Perform the computation

LRC = LRCEstimator(mu,nu,lambda,eps,verbosityFlag);

fprintf('LRC is %d \n',LRC)

