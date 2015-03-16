function [simoutput]=run_simulation_generic(parameters, simulation_options)
% function [simoutput]=run_simulation_generic(parameters, simulation_options)
% % vector of parameters corresponding to the parameter 
% % inputs in model_name:
% parameters  = [    ]; 
%
% % simulation_options 
% % data structure with the following fields:
%
%
% % name of the simulink model (e.g., filename.mdl)
% simulation_options.model_name             = 'filename';
%
%
% % This is a data structure with field names that correspond
% % to the options in simset (help simset to see the different 
% % options here). For example to specify the ode45 solver:
% simulation_options.default_simopts.Solver = 'ode45';
%
% % the default solver is ode23s 
% 
% % Simulation methodology. If you're using automatic equation generation,
% % then two methods are available for running simulations
% % 
% %   - m-file   -> Uses matlabs m-file scripting interface for running 
% %                 simulations (default method)
% % 
% %   - simulink -> uses compiled C and is faster, but you may not have a
% %                 license for it
% % This is set using: 
%
% simulation_options.integrate_with = 'm-file';
% % or
% simulation_options.integrate_with = 'simulink';
%
% % when simulating the system for generating plots it may be important
% % to force sampling at importatn output times (such as bolus events, 
% % times when the infusion rate changes, etc). Set this value to 'yes'
% % to accomplish this (default is 'no')
%
% simulation_options.include_important_output_times = 'yes';
%
% simulation_options.output_times  = 0:100;
%
% % output times of states and model outputs
% simulation_options.output_times  = 0:100;
%
% % initial conditions of the states in the model for a model with two states:
% simulation_options.initialstate = [0 0];
% 
% % Note:
% %The initial state applies to the first value of output_times. This is
% %important to consider when running the simulation for a couple reasons:
% %    o inputs (e.g., bolus values occur before measured outputs)
% %    o system isn't running at steady state and changes 
% %      before the desired output times
% % If you want an input (e.g. bolus) to occur at a time before your first 
% % measured output, it's important to make the first output time occur 
% % at the time of this first input. 
% % 
% % In some cases states like diseases will progress in the absence of
% % treatment. For example, say you are modeling anti-tumor drugs. These data
% % usually consider measurement times relative to the beginning of the study.
% % So the 'initialstate' above actually refers to the tumor volume at time
% % zero. In this example, say we have two measurements at times 10 and 20.
% % Then output_times should look something like:
% % 
% %    simulation_options.output_times  = [0 10 20];
% % 
% % This will give the correct simulated results, and the first values of the
% % output value below can be discarded
%  
% %
% % Model Inputs: 
% %     - constant infusions (infusion_rates)
% %     - bolus inputs (bolus_inputs)
% %     - constant values (myconst)
% %     - time varying/piece-wise continuious values (mytimevarying)
% 
% % Field: infusion_rates
% 
% % infusion_rates inputs that shift from one level to another (e.g., shifting 
% % from no drug to a continuous infusion of 20 mg/kg for a period of time and 
% % then back to zero. This field is an array with two fields.   
% % The 'name' field is a string containing the name the simulink model
% % references for this input. The 'levels' field contains the time values 
% % in the first row, and the second row contains the values the input
% % switches to at the corresponding times
% simulation_options.infusion_rates(1).name = 'druginput';    % name used in simulation
% simulation_options.infusion_rates(1).levels = [ 0 10 15 20  % time  values
%                                                 1 0  2  0]; % input level
% % Simulink block "From Workspace", infusion_rates.druginput
%
% % In this example the input level will step up to 1 at time zero and return to
% % zero at time 10. it will then step up to 2 at time 15 and back down to zero
% % at time 20;
%
% % Field: bolus_inputs   
%
% % Matrix describing bolus events. 
% %  - The first row represents the times the bolus occurs
% %  - The second row represents the corresponding bolus values
% %  - The thrid row represents the state the bolus value is added to
% % For bolus values in mulitple states simply add two rows for each state
% simulation_options.bolus_inputs = [ 0  1 3   5  9  13 20      % bolus times
%                                     5  0 11  3  0   4  0      % bolus values for state 1
%                                     1  1  1  1  1   1  1      % state      
%                                     0 10  0  0  2   0  3      % bolus values for state 2
%                                     2  2  2  2  2   2  2]     % state      
%
% % Field: myconst        
% simulation_options.myconst(1).name   = 'constname';    % name used in simulation
% simulation_options.myconst(1).value  = 10;             % value for that constant
%       this is accessible from from simulink by creating a constant block and
%       using 'myconst.constname' as the value
% % Simulink block "Constant", myconst.constname
%
%
%
% % Field: timevarying
% simulation_options.mytimevarying(1).name   = 'tvname';          % name used in simulation
% simulation_options.mytimevarying(1).value  = [ 0 1 2 3 4 5 6    % times
%                                                0 1 2 3 3 3 0 ]; % values at those times
% % Simulink block "From Workspace", mytimevarying.tvname
%
% 
%
%   simoutput -- data structure with the following fields:
%       t -- output times
%       x -- states at corresponding times
%       y -- matrix of outputs with the following column values:
% 
% % Note: if you needed to add an extra output time to account for bolus events
% % or nonsteady-state systems as mentioned above, this extra information can be
% % removed in the following manner:
% %  simoutput.t = simoutput.t(2:end);
% %  simoutput.x = simoutput.x(2:end,:);
% %  simoutput.y = simoutput.y(2:end,:);
%
%
%
% % Alternative input format: 
% % the variable all_inputs is created which concatenates parameters, constants,
% % infusion rates, and timevarying inputs (in that order). To use this format,
% % crate a "From Workspace" block in Simulink and make this the _only_ input
% % into the mux block. Set the dimension of the input to the mux block -1 to
% % have that automatically selected. The u[] inputs in C have to be constructed
% % in the following order 
% %  u[0]     to u[p-1]       represents the parameters         (p=> number of parameters)
% %  u[p]     to u[p+c-1]     represents the constant inputs    (c=> number of constant inputs)
% %  u[p+c]   to u[p+c+r-1]   represents the rate inputs        (r=> number of infusion rate inputs) 
% %  u[p+c+r] to u[p+c+r+v-1] represents the timevarying inputs (v=> number of time varying inputs) 
%
%
%
%  version devel
%     * Added fail safe in simulink selection. This way if simulink is
%       selected and it fails (because the toolbox isn't available) the
%       simulation will automatically be run as an m-file
%
%  version 2012.06
%     * Abstracted out the simulation running portion to allow two options:
%        - running simulation with simulink (previous)
%        - running the simulation using mscript (current)
%
%  version 2012.03
%     * Adjusted internal sampling strategy, to make better use of resources
%       
%     * Added all_inputs to eliminate the need to edit simulink, the
%       individual inputs (constant, infusion, timevarying) are still
%       available for flexability
%
%  version 2012.02
%     * Adjusted sampling leading up to bolus events
%
%

% JMH check timevarying to make sure it still works after all the different
% modifications that have been made

warning off all;

% We need to combine the output_times vector with other time vectors
% internally. This is done to ensure the ode solver evaluates at important
% points (when infusions begin) as well as the desired outputs. In order to
% maintain internal consistency, we're going to make sure the output_times is
% a column vector regardless of what was specified by the users

if(isrow(simulation_options.output_times))
   simulation_options.output_times = simulation_options.output_times';
end

% It can be important to force the solver to evaluate 
% the system at specific times to make sure all events 
% are observed. The way bolus values are handled means 
% the system will be evaluated at each bolus event. However
% other events must be accounted for explicitly. This includes 
% the time varying inputs like infusion_rates and timevarying.
% The times these events occur are stored in the 
% 'important_times' variable
%
important_times = [simulation_options.output_times'];

%
% Setting the default options
%
% ode solver
if not(isfield(simulation_options, 'default_simopts'))
  simulation_options.default_simopts.Solver = 'ode23s';
elseif not(isfield(simulation_options.default_simopts, 'Solver'))
  simulation_options.default_simopts.Solver = 'ode23s';
end

% integration tool
if not(isfield(simulation_options, 'integrate_with'))
  simulation_options.integrate_with = 'm-file';
end

% include important times in the output
if not(isfield(simulation_options, 'include_important_output_times'))
  simulation_options.include_important_output_times = 'no';
end


if isfield(simulation_options, 'bolus_inputs')
  % calcualting the number of bolus inputs and storing that 
  % value in num_bolus
  [nr_bolus, num_bolus] = size(simulation_options.bolus_inputs);
else
   nr_bolus  = 0;
   num_bolus = 0;
end


% Simulink block inputs will be populated below. Here we're
% just initializing the blocks data structure with pieces for 
% each type of input:

simulation_options.blocks.myconst        = struct;
simulation_options.blocks.infusion_rates = struct;
simulation_options.blocks.mytimevarying  = struct;
simulation_options.blocks.all_inputs     = [];

% creating constant input values
if(isfield(simulation_options, 'myconst'));
  for inptno=1:numel(simulation_options.myconst)
      constant_name    = simulation_options.myconst(inptno).name;
      constant_value   = simulation_options.myconst(inptno).value; 
      eval(sprintf('simulation_options.blocks.myconst.%s = constant_value;', constant_name));
  end
end

% converting the continuious inputs specified by the user into
% a vector usable by simulink
if(isfield(simulation_options, 'infusion_rates'));
  infusion_rate_names = [];
  for inptno=1:numel(simulation_options.infusion_rates)
      input_name    = simulation_options.infusion_rates(inptno).name;
      input_profile = make_infusion(simulation_options.infusion_rates(inptno).levels, ...
                                    simulation_options.output_times);

      infusion_rate_names = [infusion_rate_names, {input_name}];
      eval(sprintf('simulation_options.blocks.infusion_rates.%s = input_profile;', input_name));

      % adding the infusion rate step changes to the important_times 
      % vector
      important_times = [important_times input_profile(:,1)'];
  end
end


% creating piecewise continuious inputs
if(isfield(simulation_options, 'mytimevarying'));
  timevarying_input_names = [];
  for inptno=1:numel(simulation_options.mytimevarying)
      time_varying_value      = [];
      time_varying_name       = simulation_options.mytimevarying(inptno).name;

      % if the timevarying input does not last through the duration of the
      % simulation, then the last value is held constant 
      if((simulation_options.mytimevarying(inptno).value(1,:))<max(simulation_options.output_times))
          simulation_options.mytimevarying(inptno).value(:, end+1) = [max(simulation_options.output_times); simulation_options.mytimevarying(inptno).value(2, end) ];
      end
      time_varying_value      = [simulation_options.mytimevarying(inptno).value(1,:)' simulation_options.mytimevarying(inptno).value(2,:)'];
      timevarying_input_names = [timevarying_input_names , {time_varying_name}];
      eval(sprintf('simulation_options.blocks.mytimevarying.%s = time_varying_value;', time_varying_name));
      % adding the timevarying time values the important_times 
      % vector
      important_times = [important_times simulation_options.mytimevarying(inptno).value(1,:)];
  end
end

if(isfield(simulation_options, 'bolus_inputs'));
    important_times = [important_times simulation_options.bolus_inputs(1,:)];
end


% Because some of the 'important' times may be repeats
% we're pulling out just the unique ones. We're also 
% sorting those times.

if(~isempty(important_times))
  important_times = sort(unique(important_times));
end

% 
% creating all_inputs, it has the form of
%
%  [important_times parameters constant_inputs infusion_rates piecewise_continuous_inputs]
%

% setting the times
simulation_options.blocks.all_inputs(:,1) = important_times';
%
% all_inputs: adding parameter inputs
%
offset = 1;
for inptno=1:numel(parameters)
    simulation_options.blocks.all_inputs(:,inptno+offset) = parameters(inptno);
end
% shifting over to the columns after the parameters
offset = offset+length(parameters);
%
% all_inputs: adding constant inputs
%
%keyboard
if(isfield(simulation_options, 'myconst'));
  for inptno=1:numel(simulation_options.myconst)
    % because constants can be vectors, we have to add each
    % element as a column in all_inputs
    for valno = 1:length(simulation_options.myconst(inptno).value)
     simulation_options.blocks.all_inputs(:,inptno+offset) = simulation_options.myconst(inptno).value(valno);
     % shifting over to the columns after the constant inputs
     offset = offset+1;
    end
  end
end 
%
% all_inputs: adding rate inputs
%
if(isfield(simulation_options, 'infusion_rates'));
  for inptno=1:numel(infusion_rate_names)
      input_name    = char(infusion_rate_names(inptno));
      input_profile = getfield(simulation_options.blocks.infusion_rates, input_name); 
      % interpolating the infusion rate at the sample times 
      %(first caolumn of all_inputs)
      simulation_options.blocks.all_inputs(:,inptno+offset) = interp1(input_profile(:,1), input_profile(:,2), simulation_options.blocks.all_inputs(:,1), 'linear', 'extrap');
  end
  offset = offset+length(infusion_rate_names);
end

%
% all_inputs: adding time varying inputs 
%
if(isfield(simulation_options, 'mytimevarying'));
  for inptno=1:numel(timevarying_input_names)
      input_name    = char(timevarying_input_names(inptno));
      input_profile = getfield(simulation_options.blocks.mytimevarying, input_name); 
      simulation_options.blocks.all_inputs(:,inptno+offset) = interp1(input_profile(:,1),input_profile(:,2), simulation_options.blocks.all_inputs(:,1));
  end
  offset = offset+length(timevarying_input_names);
end

if(strcmp(simulation_options.include_important_output_times, 'yes'));
   simulation_options.output_times = important_times';
end


% variables to store simulation outputs
tout = [];
xout = [];
yout = [];
                                                                                                      


%    here we step through each bolus value
%    in general we perform the following
%     o  integrate to a bolus event
%     o  stop the simulation and alter 
%        the initial condition of the state 
%        where the bolus is administered
%     o  repeat until we've administered all 
%        of the bolus values



if num_bolus > 0
  % this accounts for simulations which have bolus inputs

  % if the first bolus isn't applied until after the first output time, we have 
  % to integrate up to that bolus first bolus event
  if(simulation_options.bolus_inputs(1,1) > min(simulation_options.output_times))

    initialstate = simulation_options.initialstate;
   %options          = simset(default_simopts, ...
   %                          'SrcWorkspace', 'current', ...
   %                          'InitialState', simulation_options.initialstate);

    %pulling the output times from the fist time pout to the bolus event
    output_times = linspace(min(simulation_options.output_times),  simulation_options.bolus_inputs(1,1),10)';

    % Adding a sample point the next bolus event to prevent linear
    % interpolation between the sample point just before the bolus
    output_times = [output_times(1:end-2); (output_times(end)-500*eps); output_times(end)];

    % 
    % adding important_times to output_times 
    % 
    if(~isempty(important_times))
      output_times = include_important(output_times, important_times);
    end

    % running the simulation for this chunk
    [ttmp,xtmp,ytmp] = mysim(simulation_options, output_times, initialstate);

    % appending the simulated variables
    tout                            = [tout; ttmp(1:end-1,:)];
    xout                            = [xout; xtmp(1:end-1,:)];
    yout                            = [yout; ytmp(1:end-1,:)];

    %perserving the final state as the new initial state for the bolus inputs
    %below
    simulation_options.initialstate = xtmp(end,:);
  end

  % next we iterate through the bolus events
  for bolus_id = 1:num_bolus 
     % only simulate the bolus effects if it occurs before the last time pont
     if(simulation_options.bolus_inputs(1,bolus_id) < max(simulation_options.output_times))
        % JMH DELETE value of the current bolus
        % JMH DELETE bolus_value = simulation_options.bolus_inputs(2,bolus_id);
        
        % integrating from the current bolus to the next
        if bolus_id < num_bolus
          output_times = linspace(simulation_options.bolus_inputs(1,bolus_id), simulation_options.bolus_inputs(1,bolus_id+1),10)';
          % Adding a sample point the next bolus event to prevent linear
          % interpolation between the sample point just before the bolus
          output_times = [output_times(1:end-2); (output_times(end)-500*eps); output_times(end)];
        else
          output_times = linspace(simulation_options.bolus_inputs(1,bolus_id),  max(simulation_options.output_times),10)';
        end


        
        % initializing system
        if bolus_id == 1;
          % first bolus
          initialstate = simulation_options.initialstate;
        else
          % bolus_id^th bolus
          initialstate = xtmp(end,:);
        end
        
        % adding bolus values to their respective compartments
        % this is the number of compartments with bolus information
        num_bolus_compartments = (nr_bolus - 1)/2;
        for compartment = 1:num_bolus_compartments
           tmp_bolus_compartment = simulation_options.bolus_inputs((2*compartment+1),bolus_id);
           tmp_bolus_value       = simulation_options.bolus_inputs((2*compartment),bolus_id);
        
           initialstate(tmp_bolus_compartment) = initialstate(tmp_bolus_compartment) + tmp_bolus_value;
        end
        
        % 
        % only including important_times if they exist
        % 
        if(~isempty(important_times))
          output_times = include_important(output_times, important_times);
        end
        % running the simulation for this chunk
        [ttmp,xtmp,ytmp] = mysim(simulation_options, output_times, initialstate);
        
        
        % appending the chunk of time, state, and outputs to 
        % the appropriate variables
        if bolus_id < num_bolus
           tout = [tout; ttmp(1:end-1,:)];
           xout = [xout; xtmp(1:end-1,:)];
           yout = [yout; ytmp(1:end-1,:)];
        else
           tout = [tout; ttmp(:,:)];
           xout = [xout; xtmp(:,:)];
           yout = [yout; ytmp(:,:)];
        end
     end
  end

  % Interpolating to get the outputs and states 
  % at the output times of interest
  simoutput.y = interp1(tout, yout, simulation_options.output_times, 'linear', 'extrap');
  simoutput.x = interp1(tout, xout, simulation_options.output_times, 'linear', 'extrap');
  simoutput.t = simulation_options.output_times;

 %JMH for debugging:
 %if(sum(sum(isnan(simoutput.x))) > 1)
 %    save stuff.mat;
 %end

else

  % here we have no bolus inputs so we just integrate over
  % the specified time interval assuming that there are 
  % no inputs or inputs are specified in arrays:
  %     simulation_options.infusion_rates
  %     simulation_options.myconst 
  %

  % 
  % only including important_times if they exist
  % 
  if(~isempty(important_times))
    all_times    = include_important(simulation_options.output_times, important_times);
  else
    all_times    = simulation_options.output_times;
  end

  initialstate = simulation_options.initialstate;
 %options      = simset(default_simopts, ...
 %                      'SrcWorkspace', 'current', ...
 %                      'InitialState', initialstate);
  [tout ,xout , yout ] = mysim(simulation_options,   all_times, initialstate);

  simoutput.y = interp1(tout, yout, simulation_options.output_times);
  simoutput.x = interp1(tout, xout, simulation_options.output_times);
  simoutput.t = simulation_options.output_times;
end


function [infusion]=make_infusion(levels,output_times)

switch_times  = levels(1,:);
magnitude     = levels(2,:);

delta         = 500*eps;

for i = 1:numel(switch_times)
    if (1==i)
      infusion = [switch_times(i) magnitude(i)];
    else
      infusion = [infusion
                   switch_times(i)          magnitude(i-1)
                  (switch_times(i)+delta)   magnitude(i)];
    end
end

% extending the last infusion rate to the end of 
% simulation time
if(max(output_times) > max(infusion(:,1)))
  infusion = [infusion
              max(output_times) infusion(end,2)];
end


function [simopts]=make_simopts(default_simopts);

tmp_fields  =  fieldnames(default_simopts);
simopts     =  simset();

for i=1:numel(tmp_fields);
  simopts = simset(simopts, char(tmp_fields(i)), getfield(default_simopts, char(tmp_fields(i))));
end



function [all_times] = include_important(output_times, important_times)
%
% For any given set of output_times we want to make sure that the important
% times that fall within the range of output_times are present within the
% output_times. So for example if we have the following output_times, and
% important_times
%
%  output_times    = [4.0000    4.3333    4.6667    5.0000]';
%
%  important_times = [3.0000    4.5000    5.0000    6.0000]';
%
%  all_times should return the following:
%
%  all_times  -->  [4.0000    4.3333    4.5000    4.6667    5.0000]';
%

%
% making sure important_times is a column vector
%
if(isrow(important_times))
    important_times = important_times';
end

all_times = output_times;

tmax = max(output_times);
tmin = min(output_times);

% we only want to include those times which are within the range specified in
% output_times
current_important_times = important_times(((important_times >= tmin) &  (important_times <= tmax)));

all_times = [current_important_times; output_times];

all_times = sort(unique(all_times));



%
% wrapper to switch between simulink and m-file script
%
function [t, x, y] = mysim(simulation_options,  output_times, options)



if(strcmp(simulation_options.integrate_with, 'm-file'))
  [t,x,y] = sim_m_file(simulation_options, output_times, options);

elseif(strcmp(simulation_options.integrate_with, 'simulink'))
  try
    [t,x,y] = sim_simulink(simulation_options, output_times, options);
  catch
    disp('#> Simulink failed switching to m-file');
    [t,x,y] = sim_m_file(simulation_options, output_times, options);
  end
end



%
% running the simulation with an m-file
%
function [t,x,y] = sim_m_file(simulation_options, output_times, initialstate)

[t, x, y] = auto_sim(simulation_options, output_times, initialstate);


%
% running the simulation using Simulink
%
function [t,x,y] = sim_simulink(simulation_options, output_times, initialstate)


% estracting the simulation options specified by the user
if(isfield(simulation_options, 'default_simopts'));
  default_simopts = make_simopts(simulation_options.default_simopts);
else
  default_simopts = simset; 
end


% pulling out the 
simset_options   = simset(default_simopts, ...
                      'SrcWorkspace', 'current', ...
                      'InitialState', initialstate);

%
% Defining the components that 
% are fed into simulink
%
block_names = fieldnames(simulation_options.blocks);

for block_idx = 1:length(block_names)
    eval(sprintf('%s = simulation_options.blocks.%s ;', block_names{block_idx}, block_names{block_idx}));
end

% running the simulation
[t,x,y] = sim(simulation_options.model_name, output_times, simset_options);


function [expanded_times] = expand_about_times(time_vector, output_times)
% function [expanded_times] = expand_about_times(time_vector, output_times)
%
% This function takes a vector of times and returns an expanded vector
% (expanded_times) with linearly spaced elements +/- 2 % of the total time
% range

expanded_times = [];

time_range = max(output_times) - min(output_times) ;
time_delta = time_range*.02;

for i=1:length(time_vector)
    expanded_times = [expanded_times linspace((time_vector(i)-time_delta), (time_vector(i)+time_delta), 100)];
end

