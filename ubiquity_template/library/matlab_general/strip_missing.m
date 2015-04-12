function [sample_time_stripped output_stripped]=strip_missing(sample_time, output, missing_number)
% function [sample_time_stripped output_stripped]=strip_missing(sample_time, output, missing_number)
%
% For a vector of outputs (output) and corresponding sample times
% (sample_time), missing values are normally represented by -1. This function 
% removes the missing values and corresponding times and places those in the
% corresponding variables output_stripped and sample_time_stripped,
% respectively.
%
% % example
% output      = [1  2  3  4 -1  6 -1  8]';
% sample_time = [1  2  3  4  5  6  7  8]';
%
% [sample_time_stripped output_stripped]=strip_missing(sample_time, output)
%
% % By default all missing values are removed. If you have negative values in
% % your dataset, you can use missing_number to indicate that a specific value
% % has been selected as a missing value  (e.g. -1e10)


if((exist('missing_number', 'var')))
  output_stripped      = output(output~=missing_number);
  sample_time_stripped = sample_time(output~=missing_number);
else
  output_stripped      = output(output >= 0 );
  sample_time_stripped = sample_time(output >= 0);
end


