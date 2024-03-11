function IterSettings = IterSetting( UserSettings )
% Setting iteration and convergence criterion here.
% further implementation may allow reading user-defined settings from
% file directly.

%% default value

IterSettings.autosave_period = Inf; % For what period current process will be saved automatically.
IterSettings.plot_period = 10;

IterSettings.dt = 0.0025;           % only used in transient solver
IterSettings.tf = 10;             % only used in transient solver
% inner_iter_max = 2; % only used in transient solver

end

