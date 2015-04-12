function varargout = model_gui(varargin)
% MODEL_GUI MATLAB code for model_gui.fig
%      MODEL_GUI, by itself, creates a new MODEL_GUI or raises the existing
%      singleton*.
%
%      H = MODEL_GUI returns the handle to a new MODEL_GUI or the handle to
%      the existing singleton*.
%
%      MODEL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODEL_GUI.M with the given input arguments.
%
%      MODEL_GUI('Property','Value',...) creates a new MODEL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before model_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to model_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help model_gui

% Last Modified by GUIDE v2.5 14-Nov-2012 10:23:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @model_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @model_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before model_gui is made visible.
function model_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output_graph args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to model_gui (see VARARGIN)


% setting the paths
% this only happens in the undeployed 
% case (i.e. when we're not a stand alone 
% execuitable) 
if(~isdeployed && ~ismcc)
   provision_workspace;
end


% Initializing the branding of the figure
%brand_figure(handles);

% creating the log file
OF = fopen('model_gui_log.txt', 'w');
fclose(OF);
%initializing the the log file
add_to_logfile(sprintf('%s\n', datestr(now)));
%add_to_logfile(sprintf('%s\n', pwd));

% Choose default command line output_graph for model_gui
handles.output_graph = hObject;

% Update handles structure
guidata(hObject, handles);

% populating the different elements of the gui
reset_gui(handles);
status_message(handles, 'Running simulations with initial parameter set, be patient.');
% forcing the GUI to update different things
drawnow;
manage_state(handles);
update_simulation_output_figure(handles);
status_message(handles, 'A/C PDM GUI Initialized');

function [] = update_simulation_output_figure(handles)
 
% pulling up the system state to see if it is 'stale' 
cfg = manage_state(handles);
if(strcmp('stale', cfg.options.gui_state))
  status_message(handles, 'Running simulation, be patient');
  exit_status = run_simulation(handles);
  if(strcmp(exit_status, ''))
    status_message(handles, 'Simulation complete, updating figure');
    % reloading the system state to 
    % reflect the simulation 
    cfg = manage_state(handles);
  else
    status_message(handles, exit_status);
  end
else
 status_message(handles, 'Replotting from previous simulation');
end


axes(handles.plot_output);
cla;
hold on;
markers = [{'b-'}
           {'g-'}
           {'r-'}
           {'k-'}
           {'b-.'}
           {'g-.'}
           {'r-.'}
           {'k-.'}];

possible_outputs = get(handles.output_list, 'string');
selected_outputs = get(handles.output_list, 'value');

legend_handles = [];
legend_text    = [];


% defaulting the plot_times to sim_time
plot_times = cfg.simout_mapped.times.sim_time;

% now if a timescale was specified, we pick that up
 if(isfield(cfg.options, 'misc'))
   if(isfield(cfg.options.misc, 'TS'))
     eval(sprintf('plot_times = cfg.simout_mapped.times.%s;', cfg.options.misc.TS ));
   end
 end
marker_idx = 1;
for selected_idx = 1:length(selected_outputs)
  marker = char(markers(marker_idx));
  output_name =char(possible_outputs(selected_outputs(selected_idx)));
  
  current_output = getfield(cfg.simout_mapped.outputs, output_name);
  
  % plotting each selected output
  htmp = plot(plot_times, current_output, marker);
  
  % storing legend handles and text
  legend_handles = [legend_handles; htmp];
  legend_text    = [legend_text; {output_name}];  
  
  % making sure we don't go past the number of markers available 
  marker_idx = marker_idx + 1;
  if(marker_idx > length(markers))
    marker_idx = 1;
  end
end


   %tightening the axis
  axis tight;
  %disp('here')
   % now overwriting with the prespecified values

   %if(isfield(cfg.options.misc, 'xlim'))
     set(gca, 'xlim', [0 str2num(get(handles.output_times, 'string'))]);
   %end
   if(isfield(cfg.options.misc, 'ylim'))
     set(gca, 'ylim', str2num(cfg.options.misc.ylim));
   end
legend(legend_handles, legend_text,  'Interpreter', 'none')


% checking for guide information
if(isfield(cfg.options.misc, 'guides'))
  guide_names = fieldnames(cfg.options.misc.guides);
  
  current_xlim = get(gca, 'xlim');
  current_ylim = get(gca, 'ylim');
  
  for guide_idx = 1:length(guide_names)
    current_guide = getfield(cfg.options.misc.guides, char(guide_names(guide_idx)));
    if(strcmp(char(current_guide.type), 'horizontal'))
      plot(linspace(current_xlim(1), current_xlim(2), 40), ...
           ones(1,40).*current_guide.level, ...
           current_guide.marker, ...
           'color', current_guide.color);
         
    elseif(strcmp(char(current_guide.type), 'vertical'))
      plot(ones(1,40).*current_guide.level, ...
           linspace(current_ylim(1), current_ylim(2), 40), ...
           current_guide.marker, ...
           'color', current_guide.color);
    end
    
  end
  
end

prepare_figure('present');

if(isfield(cfg.options.misc, 'TS'));
  xlabel(sprintf('time (%s)', cfg.options.misc.TS))
else
  xlabel('time');
end




function [SIMINT_failure_message] = run_simulation(SIMINT_handles)
SIMINT_failure_message = '';
SIMINT_simout_mapped = {};

% getting the current state 
% of the gui
cfg = manage_state(SIMINT_handles);


% setting up basic simulation options
cfg.options.simulation_options.model_name             = 'ode_simulation';
cfg.options.simulation_options.include_important_output_times = 'yes';


% defining parameters
% default values set in reset_gui
parameters =  cfg.parameters.values;


%-----------------------------------------------------------------------------------------
% Pulling bolus information from the GUI
%
if(bolus_events_exist(cfg))
  try
    SIMINT_MY.dosing_number        = str2num(get(SIMINT_handles.dosing_number, 'string'));
    SIMINT_MY.dosing_time          = str2num(get(SIMINT_handles.dosing_time,   'string'));
    SIMINT_MY.repeat_dose_checkbox = get(SIMINT_handles.repeat_dose_checkbox, 'value');
    try
    SIMINT_MY.tmp_data = get(SIMINT_handles.bolus_table, 'Data');
    
    % Defaulting to the 
    cfg.options.inputs.bolus.times.values =  str2num(cell2mat(SIMINT_MY.tmp_data(1,2)));
    if(SIMINT_MY.repeat_dose_checkbox == 1)
      % If repeat dosing is selected then we create the time vector with the
      % entered values plus the repeat values
      SIMINT_MY.repeat_dose_times    = [0:(SIMINT_MY.dosing_number-1)].*SIMINT_MY.dosing_time+SIMINT_MY.dosing_time+(cfg.options.inputs.bolus.times.values(1,end));
      cfg.options.inputs.bolus.times.values = [];
      cfg.options.inputs.bolus.times.values = [str2num(cell2mat(SIMINT_MY.tmp_data(1,2))), SIMINT_MY.repeat_dose_times];
    end
 
    % Subsequent rows contain dosing information
    % removing the first row from tmp_data;
    SIMINT_MY.tmp_data = SIMINT_MY.tmp_data(2:end,:);
    SIMINT_MY.bolus_states = SIMINT_MY.tmp_data(:,1);
    % for each state
    for SIMINT_state_idx = 1:length(SIMINT_MY.bolus_states)
      SIMINT_MY.offset          = (SIMINT_state_idx)*2;
      % getting the loading dose vector specified in the gui  
      SIMINT_MY.mag_vector      = str2num(cell2mat(SIMINT_MY.tmp_data(SIMINT_state_idx,2)));
      
      % adding repeating dose magnituides here for the final dose level
      if(SIMINT_MY.repeat_dose_checkbox == 1)
        SIMINT_MY.mag_vector = [SIMINT_MY.mag_vector, SIMINT_MY.mag_vector(end).*ones(1,SIMINT_MY.dosing_number)];
      end

      % Adding the magnitude vector to cfg:
      eval(sprintf('cfg.options.inputs.bolus.species.%s.values = SIMINT_MY.mag_vector;', SIMINT_MY.bolus_states{SIMINT_state_idx}));
    end
    catch error_message
      add_exception_to_logfile(error_message);
      SIMINT_failure_message = strcat(SIMINT_failure_message, '*unable to add bolus events for compartments*');
    end
    catch error_message
    add_exception_to_logfile(error_message); 
    SIMINT_failure_message = strcat(SIMINT_failure_message, '*unable to evaluate bolus times*');
  end
end
%-----------------------------------------------------------------------------------------



%-----------------------------------------------------------------------------------------
% Pull infusion information from the GUI
if(infusion_rates_exist(cfg))

  clear SIMINT_MY;
  SIMINT_MY.tmp_data = get(SIMINT_handles.infusion_table, 'Data');
  %tmp_data is a cell array with the following format:
  %  'myrate1'    'times'     '-100 0 4'    'hours'
  %  ''           'levels'    '0 1 0'       'mg/hr'
  %  'myrate2'    'times'     '-100 0 4'    'hours'
  %  ''           'levels'    '0 3 0'       'mg/hr'
  for SIMINT_MY_INFUSION_IDX = 1:(length(SIMINT_MY.tmp_data(:,1))./2)
    
    SIMINT_MY.offset = (SIMINT_MY_INFUSION_IDX-1)*2 + 1;
    SIMINT_MY.rate_name = SIMINT_MY.tmp_data{(SIMINT_MY.offset),  1};
    SIMINT_MY.times     = SIMINT_MY.tmp_data{(SIMINT_MY.offset),  3};
    SIMINT_MY.levels    = SIMINT_MY.tmp_data{(SIMINT_MY.offset+1),3};
    eval(sprintf('cfg.options.inputs.infusion_rates.%s.times.values  = [%s];', SIMINT_MY.rate_name, SIMINT_MY.times));
    eval(sprintf('cfg.options.inputs.infusion_rates.%s.levels.values = [%s];', SIMINT_MY.rate_name, SIMINT_MY.levels));
  end
end
%-----------------------------------------------------------------------------------------


%-----------------------------------------------------------------------------------------
% setting up output times 
try 
 % default no scaling of time
 SIMINT_MY.sim_timescale = 1;
 SIMINT_MY.TS            = 'sim_time';

 % the time horizon shown in the figure is scaled to observable time, so we
 % have to scale it back to sim_time
 % if there one specified then we overwrite it
 if(isfield(cfg.options, 'misc'))
   if(isfield(cfg.options.misc, 'TS'))
      SIMINT_MY.sim_timescale = getfield(cfg.options.time_scales, cfg.options.misc.TS);
   end
 end

 
 eval(sprintf('cfg.options.simulation_options.output_times = linspace(0,%s,2000)./%.5e;',cfg.options.plot.final_time, SIMINT_MY.sim_timescale));

 if(isrow(cfg.options.simulation_options.output_times))
   cfg.options.simulation_options.output_times  = cfg.options.simulation_options.output_times';
 end
 catch error_message
   add_exception_to_logfile(error_message); 
  SIMINT_failure_message = strcat(SIMINT_failure_message, '*unable to assign output times*');
end
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% checking to see if a gui callback was specified
% if it was we're going to try to execute it here
if(isfield(cfg.options, 'misc'))
  if(isfield(cfg.options.misc, 'callback_gui'))
    try
      eval(cfg.options.misc.callback_gui);
    catch error_message
     add_exception_to_logfile(error_message); 
     SIMINT_failure_message = strcat(SIMINT_failure_message, '*gui callback failed*');
    end
  end
end
%-----------------------------------------------------------------------------------------


%trying to run the simulation
try
  % simulating the system
  try
    % first we try ode15s
    SIMINT_simulation_options.default_simopts.Solver = 'ode15s';
    SIMINT_simout_mapped  = run_simulation_ubiquity(parameters, cfg);
  catch
    % if that fails then we try ode23s
    SIMINT_simulation_options.default_simopts.Solver = 'ode23s';
    SIMINT_simout_mapped  = run_simulation_ubiquity(parameters, cfg);
  end
  manage_state(SIMINT_handles, SIMINT_simout_mapped);
catch error_message
  add_exception_to_logfile(error_message); 
  SIMINT_failure_message = strcat(SIMINT_failure_message, '*execution of run_simulation_ubiquity failed*');
end




function [] = reset_gui(handles, selected_parameter_set_idx)
% clearing any previous gui_state
% this forces the get_current_state call
% to create a default state
if(exist('gui_state.mat', 'file'))
  delete gui_state.mat;
end

% if no parameter set is defined, we define it as 1 
% (the default parameter set)
if(~exist('selected_parameter_set_idx', 'var'))
  selected_parameter_set_idx = 1;
end

cfg = manage_state(handles, [], selected_parameter_set_idx);

%
% Populating the pulldown menu for the different systems
% 

psdata = cell(0);
for set_idx = 1:length(cfg.options.mi.parameter_sets_reverse)
  parameter_set_id   = cfg.options.mi.parameter_sets_reverse{set_idx};
  parameter_set_name = getfield(getfield(cfg.parameters.sets, parameter_set_id), 'name');
  psdata(set_idx)    = {parameter_set_name};
end
set(handles.parameter_set_popupmenu, 'String', psdata);

% Making sure the correct parameter set is selected in the GUI:
set(handles.parameter_set_popupmenu, 'Value', selected_parameter_set_idx);

%
% parameters
%
ptypes = unique(cfg.parameters.type);
pdata = [];
for type_idx = 1:length(ptypes)
  % this is the current parameter type
  tmp_type = char(ptypes(type_idx));  
  found_type = false;
  for parameter_idx = 1:length(cfg.parameters.names)
     tmp_pname = char(cfg.parameters.names(parameter_idx));
     % checking to see if the current parameter is in the correct type
     % and seeing if it is editable
     if(strcmp(tmp_type, char(cfg.parameters.type(parameter_idx))) & ...
        strcmp(char(cfg.parameters.editable(parameter_idx)), 'yes'))
       if(not(found_type))
         %Creating a heading row
         pdata      = [pdata; {tmp_type}, {''}, {''}];
         found_type = true;
       end
       % Only adding the parameter if it has the same type as
       % the current type grouping specified in tmp_type
       pdata = [pdata; {tmp_pname} {mynum2str(cfg.parameters.values(parameter_idx))} cfg.parameters.units(parameter_idx)];
     end
   end
end
set(handles.parameters_table, 'Data', pdata);
set(handles.parameters_table, 'ColumnName', [{'Parameter'}, {'Value'}, {'Units'}]);
set(handles.parameters_table, 'ColumnEditable', logical([0 1 0]));

%
% resetting the different input fields (bolus and infusion rates)
%
set(handles.bolus_table, 'ColumnName', [{''}, {'Value'}, {'Units'}]);
if(bolus_events_exist(cfg))
  % getting the compartments that can receive bolus values
  compartment_names = fieldnames(cfg.options.inputs.bolus.species);
  tmp_data = [{'dose times'}, {num2str(cfg.options.inputs.bolus.times.values)} {cfg.options.inputs.bolus.times.units}];
  for compartment_idx = 1:length(compartment_names)
    compartment = char(compartment_names(compartment_idx));
    dosing = getfield(cfg.options.inputs.bolus.species, compartment);
    tmp_data = [tmp_data; (compartment), {num2str(dosing.values)} {dosing.units}];
    
  end
  set(handles.bolus_table, 'Data', tmp_data);
  set(handles.bolus_table, 'ColumnEditable', logical([0 1 0]));
  try
    set(handles.dosing_time_units, 'String', cfg.options.inputs.bolus.times.units)
  end
else
  % disabling the bolus fields
  set(handles.bolus_table,   'enable', 'off');
  set(handles.dosing_time,   'enable', 'off');
  set(handles.dosing_number, 'enable', 'off');
end

% if there is no model diagram file (system.jpg, or system.png) 
% we disable the 'model structure' button
if(~exist('system.jpg', 'file') & ~exist('system.png', 'file'))
  set(handles.view_model_button, 'enable', 'off');
end


set(handles.infusion_table, 'ColumnName', [{'Input'}, {'Components'}, {'Value'},{'Units'}]);
if(infusion_rates_exist(cfg))
  rate_names = fieldnames(cfg.options.inputs.infusion_rates);
  tmp_data = [];
  for rate_idx = 1:length(rate_names);
    rate = char(rate_names(rate_idx));
    rate_info = getfield(cfg.options.inputs.infusion_rates, rate);
    tmp_data = [tmp_data; {rate} {'times'}  {num2str(rate_info.times.values)}  {rate_info.times.units}];
    tmp_data = [tmp_data; {''}   {'levels'} {num2str(rate_info.levels.values)} {rate_info.levels.units}];
  end

  set(handles.infusion_table, 'Data', tmp_data);
  set(handles.infusion_table, 'ColumnEditable', logical([0 0 1 0]));
else
  % disabling the infusion table
   set(handles.infusion_table, 'enable', 'off');
end




% Plot output_graph fields
output_names_visible = {};
output_names         = fieldnames(cfg.options.mi.outputs);

for(output_idx = 1:length(output_names))
  if(regexp(output_names{output_idx}, '^QC'))
  else
    output_names_visible = [output_names_visible , output_names{output_idx}]; 
  end  
end
set(handles.output_list, 'string', output_names_visible);

%
% defining misc information
%

% output_graph times
 if(isfield(cfg.options, 'misc'))
   if(isfield(cfg.options.misc, 'output_times'))
     eval(sprintf('output_times = %s; ', cfg.options.misc.output_times));
     final_time = max(output_times);
     % scaling the final time to the unit specified in TS
     if(isfield(cfg.options.misc, 'TS'))
       final_time = final_time.*getfield(cfg.options.time_scales, cfg.options.misc.TS);
     end
     set(handles.output_times, 'string', num2str(final_time));
   else
     set(handles.output_times, 'string', '100');
   end
 end
 

 
 axes(handles.plot_output)
 %brand_figure(handles);
 

 
% check to see if bolus events exist in cfg
function [ic_status] = initial_conditions_exist(cfg)
ic_status = false;
if(isfield(cfg.options, 'initial_conditions'))
  ic_status = true;
end
  
  
function [bolus_status] = bolus_events_exist(cfg)
% defaulting to no bolus
bolus_status = false;
if(isfield(cfg.options, 'inputs'))
  if(isfield(cfg.options.inputs, 'bolus'))
    bolus_status = true;
  end
end
 
function [input_status] = time_varying_inputs_exist(cfg)
  input_status = false;
if(isfield(cfg.options, 'inputs'))
  if(isfield(cfg.options.inputs, 'time_varying_inputs'))
    input_status = true;
  end
end

function [infusion_status] = infusion_rates_exist(cfg)
% defaulting to no bolus
infusion_status = false;
if(isfield(cfg.options, 'inputs'))
  if(isfield(cfg.options.inputs, 'infusion_rates'))
    infusion_status = true;
  end
end

function [cfg] = fetch_previous_state(handles)
failure_message = '';
if(exist('gui_state.mat', 'file'))
   load gui_state.mat cfg;
else
  failure_message = '*unable to load previous state, returning current*';
  cfg = manage_state(handles);
end
% this function loads the previous state, 
% updates those values with the current state
% and saves that state
function [cfg] = manage_state(handles, simout_mapped, selected_parameter_set_idx)

 failure_message = '';
 if(exist('gui_state.mat', 'file'))
  load gui_state.mat cfg;
  % now we update the state using the current values in the GUI
  % pulling info from the gui
  % bolus information


  % output times
  try
    cfg.options.plot.final_time = get(handles.output_times, 'string');
   catch error_message
   add_exception_to_logfile(error_message); 
    failure_message = strcat(failure_message, '*Unable to process output times*');
  end

  % processing parameter values
  try
  pdata   = get(handles.parameters_table, 'Data');
  pnames  = pdata(:,1);
  pvalues = pdata(:,2);
  ptypes  = unique(cfg.parameters.type);
  for parameter_idx = 1:length(pnames)
    pname = char(pnames(parameter_idx));
    % checking to see if the 'name' is really a parameter type
    if(not(ismember(pname, ptypes)))
      cfg.parameters.values(getfield(cfg.options.mi.parameters, pname)) = str2num(char(pvalues(parameter_idx)));
    end
  end
   catch error_message
   add_exception_to_logfile(error_message); 
    failure_message = strcat(failure_message, '*Unable to process parameters*');
  end
 else
  % this should happen when the gui is initialized or 
  % when the reset button has been pushed
  cfg = auto_fetch_system_information;
  % selecting the defualt parameter set
  % this is overwritten in the next conditional
  % statement if a parameter set has been selected
  cfg = select_set(cfg, 'default');
  
  if(exist('selected_parameter_set_idx', 'var'))
    % pulling out selected parameter set short name:
    active_parameter_set  = cfg.options.mi.parameter_sets_reverse{selected_parameter_set_idx};
    cfg = select_set(cfg, active_parameter_set);
    % values stores the current editable and non-editable parameters
    % by default it contains the 'default' parmeter set, here we overwrite it
    % with whatever has been selected in the gui.
    cfg.parameters.values = getfield(getfield(cfg.parameters.sets, active_parameter_set), 'values');
  end
  % setting the state to stale to force an update of the
  % simulation results
  cfg.options.gui_state = 'stale';
 end
 
 if(not(strcmp(failure_message, '')))
   status_message(handles, failure_message);
 end
 
 % storing simulation output
 % Checking to see if simout_mapped exists
 if(exist('simout_mapped', 'var'))
   % now making sure it's not empty
   if(~isempty(simout_mapped))
     % this stores the simulation output
     cfg.simout_mapped     = simout_mapped;
     % this sets the gui state to 'fresh'
     % indicating that it represents the simulation
     % corresponds to the entries currently in the GUI
     % any changes to these entries (parameters, dosing, etc.) 
     % wll change this to 'stale'
     cfg.options.gui_state = 'fresh';
   end
 end
 save gui_state.mat cfg;
 
 function [] = status_message(handles, message)
   
 %current_message = get(handles.message, 'string');
 %new_message =strcat(current_message,message);
 set(handles.message, 'string', message);

 %{
 function [] = brand_figure(handles)
axes(handles.logo);
set(gca, 'ytick', [], 'xtick', [], 'xticklabel', [], 'yticklabel', []);
box on;
[imgLogo, map, imgAlpha] = imread('Logo-square-transparent.png');
[imgLogo] = imread('Logo-square-white.png');
axesChildHandle = imshow(imgLogo, 'InitialMagnification', 'fit');
hlogo = imshow(imgLogo);
set(hlogo, 'AlphaData', imgAlpha);

axes(handles.logo_full);
set(gca, 'ytick', [], 'xtick', [], 'xticklabel', [], 'yticklabel', []);
box on;
[imgLogo_full, map, imgAlpha] = imread('Logo-rectangle.png');
axesChildHandle = imshow(imgLogo_full, 'InitialMagnification', 'fit');
hlogo = imshow(imgLogo_full);
set(hlogo, 'AlphaData', imgAlpha);
%}




function [] = update_bolus_fields(handles,cfg)


% defining the bolus times and time scale
set(handles.dose_times, 'string', mat2str(cfg.options.inputs.bolus.times.values))
if(isfield(cfg.options.inputs.bolus.times, 'scale'))
  set(handles.time_scale, 'string', cfg.options.inputs.bolus.times.scale);
end 

compartments = get(handles.compartment_menu, 'string');
selected     = get(handles.compartment_menu, 'value');

compartment = char(compartments(selected));

dosing = getfield(cfg.options.inputs.bolus.species, compartment);

%pulling the dose values and scaling parameters for the compartment
set(handles.bolus_values, 'string', mat2str(dosing.values));
set(handles.bolus_scale,  'string', dosing.scale);












% --- Outputs from this function are returned to the command line.
function varargout = model_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output_graph args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output_graph from handles structure
varargout{1} = handles.output_graph;



function dose_times_Callback(hObject, eventdata, handles)
% hObject    handle to dose_times (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dose_times as text
%        str2double(get(hObject,'String')) returns contents of dose_times as a double
manage_state(handles);


% --- Executes during object creation, after setting all properties.
function dose_times_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dose_times (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%
%
% --- Executes on selection change in compartment_menu.
function compartment_menu_Callback(hObject, eventdata, handles)
% hObject    handle to compartment_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns compartment_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from compartment_menu
cfg=fetch_previous_state(handles);
update_bolus_fields(handles,cfg);
manage_state(handles);



% --- Executes during object creation, after setting all properties.
function compartment_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to compartment_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in output_list.
function output_list_Callback(hObject, eventdata, handles)
% hObject    handle to output_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns output_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from output_list



% --- Executes during object creation, after setting all properties.
function output_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bolus_scale_old_Callback(hObject, eventdata, handles)
% hObject    handle to bolus_scale_old (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bolus_scale_old as text
%        str2double(get(hObject,'String')) returns contents of bolus_scale_old as a double


% --- Executes during object creation, after setting all properties.
function bolus_scale_old_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bolus_scale_old (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







% --- Executes during object creation, after setting all properties.
function output_times_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_times (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in toggle_yscale.
function toggle_yscale_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_yscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



axes(handles.plot_output);
current_yscale = get(handles.plot_output, 'yscale');
prepare_figure('present');
if(regexpi(current_yscale, 'linear'))
  set(handles.plot_output, 'yscale', 'log');
else
   set(handles.plot_output, 'yscale', 'linear');
end
% forcing an update of the branding
%brand_figure(handles);


% --- Executes on button press in toggle_grid.
function toggle_grid_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


axes(handles.plot_output);
prepare_figure('present');
grid;
% forcing an update of the branding
%brand_figure(handles);


% --- Executes on button press in update_simulation.
function update_simulation_Callback(hObject, eventdata, handles)
% hObject    handle to update_simulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% forcing the GUI to update different things
status_message(handles, 'Updating simulation, be patient.');
drawnow;
update_simulation_output_figure(handles);


% --- Executes on button press in reset_gui_button.
function reset_gui_button_Callback(hObject, eventdata, handles)
% hObject    handle to reset_gui_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
status_message(handles, 'GUI reset called');
reset_gui(handles);

% --- Executes on button press in view_model_button.
function view_model_button_Callback(hObject, eventdata, handles)
% hObject    handle to view_model_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure(1);
% removing space around the figure
set(gcf,'defaultaxesposition',[0 0 1 1]);
% getting rid of the menubar
set(gcf, 'MenuBar', 'none');
if(exist('system.jpg', 'file'))
  model_diagram = 'system.jpg';
elseif(exist('system.png', 'file'))
    model_diagram = 'system.png';
end
[model_structure] = imread(model_diagram);
image(model_structure);
axis off;





function infusion_rates_Callback(hObject, eventdata, handles)
% hObject    handle to infusion_rates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of infusion_rates as text
%        str2double(get(hObject,'String')) returns contents of infusion_rates as a double


% --- Executes during object creation, after setting all properties.
function infusion_rates_CreateFcn(hObject, eventdata, handles)
% hObject    handle to infusion_rates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in infusion_menu.
function infusion_menu_Callback(hObject, eventdata, handles)
% hObject    handle to infusion_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns infusion_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from infusion_menu


% --- Executes during object creation, after setting all properties.
function infusion_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to infusion_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







% --- Executes during object creation, after setting all properties.
function dosing_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dosing_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







% --- Executes during object creation, after setting all properties.
function dosing_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dosing_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [number_string] = mynum2str(mynumber);

if((mynumber < .1) | ...
   (mynumber > 1000))
   number_string = sprintf('%.3e',mynumber);
 else 
   number_string = sprintf('%.3f',mynumber);
 end


function []=dump_csv(handles)
cfg = manage_state(handles);


status_message(handles, 'Writing data to CSV file, be patient');

file_base_name = sprintf('simulation_data-%s', datestr(now, 'yyyy-mm-dd-HH_MM'));


%
% Dumping the simulation information
%
all_data = [];

% inserting simtime
%all_data = make_column('sim_time', cfg.simout_mapped.times.sim_time);

% now inserting time vectors for each timescale
tmp_headers =fieldnames(cfg.simout_mapped.times) ;
for tmp_idx=1:length(tmp_headers)
  tmp_name = char(tmp_headers(tmp_idx));
  tmp_data = getfield(cfg.simout_mapped.times, tmp_name);

  all_data = [all_data  make_column(tmp_name, tmp_data)];
end

% now inserting columns for each output
tmp_headers =fieldnames(cfg.simout_mapped.states) ;
for tmp_idx=1:length(tmp_headers)
  tmp_name = char(tmp_headers(tmp_idx));
  tmp_data = getfield(cfg.simout_mapped.states, tmp_name);

  all_data = [all_data  make_column(tmp_name, tmp_data)];
end

% now inserting columns for each state
tmp_headers =fieldnames(cfg.simout_mapped.outputs) ;
for tmp_idx=1:length(tmp_headers)
  tmp_name = char(tmp_headers(tmp_idx));
  tmp_data = getfield(cfg.simout_mapped.outputs, tmp_name);

  all_data = [all_data  make_column(tmp_name, tmp_data)];
end

csv_file_name = sprintf('%s.csv', file_base_name);
cell2csv(all_data, csv_file_name);


%
% Dumping the information from the GUI
%
gui_state_file_name = sprintf('%s.txt', file_base_name);
gui_state_contents = '';

% Dumping Parameter values
 try
   pdata   = get(handles.parameters_table, 'Data');
   pnames  = pdata(:,1);
   pvalues = pdata(:,2);
   punits  = pdata(:,3);
   if(length(pnames)>0)
     gui_state_contents = sprintf('%s\nParameters\n', gui_state_contents);
     for parameter_idx = 1:length(pnames)
       if(strcmp(pvalues{parameter_idx}, ''))
         gui_state_contents = sprintf('%s%s:\n', gui_state_contents, ...
                                                          pnames{parameter_idx});
       else
         gui_state_contents = sprintf('%s%s = %s (%s)\n', gui_state_contents, ...
                                                            pad_end_string(pnames{parameter_idx}, 13), ...
                                                            pad_end_string(pvalues{parameter_idx}, 15),...
                                                            punits{parameter_idx});
       end
     end
     gui_state_contents = sprintf('%s-------------------\n', gui_state_contents);
   else
     gui_state_contents = sprintf('%s\nNo parameters found', gui_state_contents);
   end
 catch
   gui_state_contents = sprintf('%sUnable to write parameters', gui_state_contents);
 end
 
 % Dumping Bolus Information
 try
   bdata   = get(handles.bolus_table, 'Data');
   bnames  = bdata(:,1);
   bvalues = bdata(:,2);
   bunits  = bdata(:,3);
   if(length(bnames)>0)
   gui_state_contents = sprintf('%s\nBolus information\n', gui_state_contents);
     for bolus_idx = 1:length(bnames)
       gui_state_contents = sprintf('%s%s = %s  (%s)\n', gui_state_contents, ...
                                                    pad_end_string(bnames{bolus_idx}, 10), ...
                                                    pad_end_string(bvalues{bolus_idx}, 10), ...
                                                    bunits{bolus_idx});                                                 
     end
     gui_state_contents = sprintf('%s-------------------\n', gui_state_contents);

   else
     gui_state_contents = sprintf('%s\nNo bolus information was found', gui_state_contents);
   end
 catch
     gui_state_contents = sprintf('%sUnable to write bolus information', gui_state_contents);
 end
 
 %Dumping Infusion Rate Information
 
  try
   rdata   = get(handles.infusion_table, 'Data');
   rnames  = rdata(:,1);
   rcomponents = rdata(:,2);
   rvalues = rdata(:,3);
   runits  = rdata(:,4);
   if(length(rnames)>0)
   gui_state_contents = sprintf('%s\nInfusion Rate information\n', gui_state_contents);
   for rate_idx = 1:length(rnames)
     if(~strcmp(rnames{rate_idx}, ''))
        gui_state_contents = sprintf('%s%s:\n', gui_state_contents, ...
                                          rnames{rate_idx});
     end
     
     gui_state_contents = sprintf('%s%s = %s  (%s)\n', gui_state_contents, ...
                                                  pad_end_string(rcomponents{rate_idx}, 10), ...
                                                  pad_end_string(rvalues{rate_idx}, 10), ...
                                                  runits{rate_idx});   
   end
   gui_state_contents = sprintf('%s-------------------\n', gui_state_contents);

   else
     gui_state_contents = sprintf('%s\nNo Infusion Rate information was found', gui_state_contents);
   end
 catch
   gui_state_contents = sprintf('%sUnable to write Infusion Rate information', gui_state_contents);
end
 
 OF = fopen(gui_state_file_name, 'w');
 fprintf(OF, '%s', gui_state_contents);
 fclose(OF);
 

 set(gcf,'PaperPositionMode', 'auto');
 image_file_name = sprintf('%s.jpg', file_base_name);
 print(gcf, '-djpeg', image_file_name);


 
 mymessage = '';
 mymessage = sprintf('%sData, GUI state, and Figure written to files with prefix:\n%s', mymessage, file_base_name);

 

status_message(handles, sprintf('%s', mymessage));


function [mycol] = make_column(title, data);

mycol = [{title}];
for data_idx = 1:length(data)
    mycol = [mycol; data(data_idx)];
end





%JMH clutter from here to the bottom

%function infusion_times_Callback(hObject, eventdata, handles)
%% hObject    handle to infusion_times (see GCBO)
%% eventdata  reserved - to be defined in a future version of MATLAB
%% handles    structure with handles and user data (see GUIDATA)
%
%% Hints: get(hObject,'String') returns contents of infusion_times as text
%%        str2double(get(hObject,'String')) returns contents of infusion_times as a double


%% --- Executes during object creation, after setting all properties.
%function infusion_times_CreateFcn(hObject, eventdata, handles)
%% hObject    handle to infusion_times (see GCBO)
%% eventdata  reserved - to be defined in a future version of MATLAB
%% handles    empty - handles not created until after all CreateFcns called
%
%% Hint: edit controls usually have a white background on Windows.
%%       See ISPC and COMPUTER.
%if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%    set(hObject,'BackgroundColor','white');
%end
%


% --- Executes on button press in csv_button.
function csv_button_Callback(hObject, eventdata, handles)
% hObject    handle to csv_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dump_csv(handles);

%function bolus_values_Callback(hObject, eventdata, handles)
%% hObject    handle to bolus_values (see GCBO)
%% eventdata  reserved - to be defined in a future version of MATLAB
%% handles    structure with handles and user data (see GUIDATA)
%
%% Hints: get(hObject,'String') returns contents of bolus_values as text
%%        str2double(get(hObject,'String')) returns contents of bolus_values as a double
%cfg = manage_state(handles);
%update_bolus_fields(handles,cfg);
%
%% --- Executes during object creation, after setting all properties.
%function bolus_values_CreateFcn(hObject, eventdata, handles)
%% hObject    handle to bolus_values (see GCBO)
%% eventdata  reserved - to be defined in a future version of MATLAB
%% handles    empty - handles not created until after all CreateFcns called
%
%% Hint: edit controls usually have a white background on Windows.
%%       See ISPC and COMPUTER.
%if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%    set(hObject,'BackgroundColor','white');
%end




function add_exception_to_logfile(exception)
add_to_logfile(sprintf('---------------------------- \n'));
add_to_logfile(sprintf('Exception encountered: \n'));
add_to_logfile(exception.message);
for idx = 1:length(exception.stack)
    tmpstr = sprintf('%s (%d) \n', exception.stack(idx).name, exception.stack(idx).line);
    add_to_logfile(tmpstr);
end
add_to_logfile(sprintf('---------------------------- \n'));
    

function add_to_logfile(string)
OF = fopen('pdm_gui_log.txt', 'a+');
disp(string)
fprintf(OF, '%s', string);
fclose(OF);


function make_state_stale()
 load gui_state.mat cfg;
 cfg.options.gui_state = 'stale';
 %disp('made statle');
 save gui_state.mat cfg;

 function str = pad_end_string(str, maxlength)
%  function str = padstring(str, maxlength)
%
%  adds spaces to the end off the string 'str' until it is lenth
%  'maxlength'
%

if(length(str) < maxlength)
  pad_length = maxlength - length(str);
  pad_string = repmat(sprintf(' '), 1,pad_length);
  str = sprintf('%s%s', str, pad_string);
end


% --- Executes when entered data in editable cell(s) in parameters_table.
function parameters_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to parameters_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
make_state_stale();


% --- Executes when entered data in editable cell(s) in bolus_table.
function bolus_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to bolus_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
make_state_stale();

% --- Executes on button press in repeat_dose_checkbox.
function repeat_dose_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to repeat_dose_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of repeat_dose_checkbox
make_state_stale();

function dosing_time_Callback(hObject, eventdata, handles)
% hObject    handle to dosing_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dosing_time as text
%        str2double(get(hObject,'String')) returns contents of dosing_time as a double

make_state_stale();


function dosing_number_Callback(hObject, eventdata, handles)
% hObject    handle to dosing_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dosing_number as text
%        str2double(get(hObject,'String')) returns contents of dosing_number as a double

make_state_stale();


% --- Executes when entered data in editable cell(s) in infusion_table.
function infusion_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to infusion_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

make_state_stale();

function output_times_Callback(hObject, eventdata, handles)
% hObject    handle to output_times (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output_times as text
%        str2double(get(hObject,'String')) returns contents of output_times as a double

make_state_stale();


% --- Executes on selection change in parameter_set_popupmenu.
function parameter_set_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to parameter_set_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns parameter_set_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from parameter_set_popupmenu
selected_parameter_set_idx =  get(handles.parameter_set_popupmenu, 'Value');
reset_gui(handles, selected_parameter_set_idx);
manage_state(handles);
update_simulation_output_figure(handles);
make_state_stale();



% --- Executes during object creation, after setting all properties.
function parameter_set_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameter_set_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
