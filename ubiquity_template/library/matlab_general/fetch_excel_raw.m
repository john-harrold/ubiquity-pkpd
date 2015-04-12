function [raw]=fetch_excel_raw(filename, sheetname)
% function [raw]=fetch_excel_raw(filename, sheetname)
%   
%  filename = string containing the excel file
%  sheetname = name of shee to fetch
%
%  raw = cell array of raw information

[a, b, raw] = xlsread(filename, sheetname);
