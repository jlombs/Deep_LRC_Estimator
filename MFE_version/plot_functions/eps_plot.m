function eps_plot(dim,trials,steps_per_phase)
%example usage: eps_plot(10,10,1e4);

parent_dir = cd(cd('..'));
addpath(parent_dir);

body = 'cube';
% dim = 20;
actual_vol = 2^dim;
% trials = 1000;
% steps_per_phase = 5e4;

[P] = makeBody(body,dim);
vols = zeros(trials, steps_per_phase);
for i=1:trials
    [vols(i,:)] = ThreadsTest(P,[],.1,'-num_t 1 -ratio 0.6', steps_per_phase, 1);
end

errors = abs(vols - actual_vol)/actual_vol;

avg_error = mean(errors);

close all;

figure(1);
plot(avg_error);

file_string = strcat('plots/error_plot_',int2str(dim),'_',int2str(trials),'_',int2str(steps_per_phase));

print(1, file_string,'-dpdf');

save_file = strcat('saved_data/error_plot_',int2str(dim),'_',int2str(trials),'_',int2str(steps_per_phase));

save(save_file,'errors');

end