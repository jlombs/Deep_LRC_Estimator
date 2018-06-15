%wrapper for ThreadsTest to do 15th Birkhoff polytope

parent_dir = cd(cd('..'));
addpath(parent_dir);

clear;
[P]=makeBody('birkhoff',5);
[volumes]=ThreadsTest(P,[],.2,'',200,100);
save(strcat('5birkhoff',int2str(randi(1e6,1))), 'volumes');
% n = 25;
% [T,~] = qr(rand(n));
% [P]=makeBody('cube',n);
% P.A = P.A*T;
% 
% [volumes]=ThreadsTest(P,[],.2,'-walk char',500,1000);
% save(strcat(strcat('char_rotcube',int2str(n)),int2str(randi(1e6,1))),'volumes');
% fprintf('finished coord\n');
% 
% [volumes]=ThreadsTest(P,[],.2,'',500,1000);
% save(strcat(strcat('har_rotcube',int2str(n)),int2str(randi(1e6,1))),'volumes');