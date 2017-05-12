function [recs] = nm_select_records(data, records, filter)
% function [recs] = nm_select_records(data, records, filter)
%
% Takes a data set created using nm_read_data with all or a subset of the data
% (derived from data.values) in records,  and filters the data according to
% the information specified in filter and returns that match. 
%
% This only works on columns with discrete values (0, 1, 2), and does not work
% on floating point values.
%
% Multiple filters  are joined by a boolean AND. While multiple options for a
% given filter are combined using an OR.
%
% For example to extract only the observations (EVID=0) in cohort 2 of studies 1-3,
% we create the filter in the following way:
%  
%  
% myfilter.evid    = 0;     
% myfilter.cohorts = 2;
% myfilter.studies = [1:3];
% 
% This creates the following boolean relationship:
% (evid = 0) and (chhorts = 2) and ((studies = 1) or (studies = 2) or (studies = 3))
%
%
% [obs] = nm_select_records(data, data.values, myfilter);
%
% This should display the first five rows of obs
% with the appropriate header on top
% [data.header; obs(1:5, :)]
%
%



% by default all of the data is returned
recs = records;


% now we process each field in filter
% to pull out the values that match
cols = fieldnames(filter); 

for cols_idx = 1:length(cols)

   column_name = cols{cols_idx};

   if(isfield(data.column.numbers, column_name))

     % this is the column number associated with column_name above
     column_number = getfield(data.column.numbers, column_name);

     % getting the values to compare
     values = getfield(filter, column_name);

     my_cond = '';
     % we're constructing the boolean values based on the info specified in
     % the filter:
     %
     %  if the filter looked like: 
     %  myfilter.subjects = [1:3], and subjects was column 2 in the dataset
     %  then the boolean equation would look something like:
     %
     %eq([recs{:, 2}]', 1) | eq([recs{:, 2}]', 2) | eq([recs{:, 2}]', 3)
     for value_idx = 1:length(values)
       if(strcmp(my_cond, ''))
         my_cond = sprintf('eq([recs{:, %d}]'', %d)', column_number, values(value_idx));
       else
         my_cond = sprintf('%s | eq([recs{:, %d}]'', %d)', my_cond, column_number, values(value_idx));
       end
     end

     
     if(strcmp(my_cond, ''))
       disp(sprintf('> No filter fields found'));
       dbstack
     else
       % now we pull out a vector of rows which match the
       % boolean condition constructed above:
       eval(sprintf('rows  = (%s);', my_cond));
       % now reducing recs down to the rows that 
       % satisfied the condition statement
       recs = recs(rows, :);
     end

   else
     disp(sprintf('> fieldname: %s not found ignoring this entry', column_name));
     dbstack
   end


end
