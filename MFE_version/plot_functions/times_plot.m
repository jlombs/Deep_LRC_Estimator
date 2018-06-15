function times_plot(dim,trials)
%example usage: times_plot(25,3);

parent_dir = cd(cd('..'));
addpath(parent_dir);

body = 'cube';
dims = 5:5:dim;
eps1 = 0.20;
eps2 = 0.10;
times1 = zeros(length(dims),trials);
samples1 = zeros(length(dims),trials);
volumes1 = zeros(length(dims),trials);
times2 = zeros(length(dims),trials);
samples2 = zeros(length(dims),trials);
volumes2 = zeros(length(dims),trials);


for i = 1:length(dims)
   [P,p] = makeBody(body,dims(i));
   for j=1:trials
   tic;
   
   [volumes1(i,j), ~, samples1(i,j)] = Volume(P,[],eps1);
   times1(i,j) = toc;
   tic;
   [volumes2(i,j), ~, samples2(i,j)] = Volume(P,[],eps2);
   times2(i,j) = toc;
   end
   fprintf('Finished dim %d/%d\n', dims(i), dims(end));
end
close all;

figure(1);
hold on;
plot(dims,mean(samples2,2).^.5,'r');
plot(dims,mean(samples1,2).^.5);
title('Dimension vs. Hit-and-run Steps (cube)');
xlabel('Dimension');
ylabel('sqrt(# Steps)');
legend('eps=0.10', 'eps=0.20');

figure(2);
hold on;
plot(dims,mean(times2,2),'r');
plot(dims,mean(times1,2));
title('Dimension vs. Time (cube)');
xlabel('Dimension');
ylabel('Time (seconds)');
legend('eps=0.10', 'eps=0.20');

suffix = strcat(int2str(dim),'_',int2str(trials));

print(1, strcat('plots/Dimension_v_Steps_Cube', suffix),'-dpdf');
print(2, strcat('plots/Dimension_v_Time_Cube', suffix),'-dpdf');


savefile = strcat('saved_data/TimesPlotData_Dim', suffix);
save(savefile, 'times1','samples1','volumes1','times2','samples2','volumes2');

end