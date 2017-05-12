function [d] = system_fetch_data(cfg, dataset, filter, column)
% function [data] = system_fetch_data(cfg, dataset, filter, column)
%  
%  dataset - Name of dataset assigned with system_load_dataset
%  filter  - Used to identify rows to return. 
%            See help nm_select_records for details
%  column  - Column name of data to return
%


% parsing user input and pulling out the dataset
isgood = true;
if(isfield(cfg.data, dataset))
  ds = getfield(cfg.data, dataset);
  if(~isfield(ds.column.names, column))
    errormsg = sprintf('Column >%s< not found in dataset >%s<', column, dataset);
    isgood = false;
  end
else
    isgood = false;
    errormsg = sprintf('Dataset >%s< not found', dataset) ;
end


if(isgood)
  if(isempty(filter))
    values = ds.values;
  else
    values = nm_select_records(ds, ds.values , filter);
  end
  d = nm_fc(ds, values, column);
else
  vp(cfg, sprintf('------------------------------------')) 
  vp(cfg, sprintf('system_fetch_data()                 ')) 
  vp(cfg, errormsg )
  vp(cfg, sprintf('------------------------------------')) 
end

