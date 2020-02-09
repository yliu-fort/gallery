clear; clc; close all;

%% Multi-grid integrator
count = dummy_integrator(1,6,0);

function count = dummy_integrator(level,num_level,count)
if(level < 1) ,return;end
if(level > num_level) ,return;end

% Find current level by counting leading zeros
fprintf("Enter level %d\n",level);
fprintf("Apply collision on level %d\n",level);
if(level < num_level)
    fprintf("Apply F-redistribution on level %d\n",level);
end

% launch new functions
count = count + 1;
count = dummy_integrator(level+1,num_level,count);
count = dummy_integrator(level+1,num_level,count);

fprintf("Apply streaming on level %d\n",level);
if(level < num_level)
    fprintf("Apply C-redistribution on level %d\n",level);
end
fprintf("Leave level %d\n",level);

end