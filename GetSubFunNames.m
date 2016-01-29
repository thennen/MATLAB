function foo = GetSubFunNames(FIDSlurp,CONTENTS)
% GETSUBFUNNAMES Get the names of the file's
% subfunctions
%
% USAGE:
%   FID = fopen('xlsread.m','rt');
%   a = GetSubFunNames(FID)

TRUE  = logical(1);
FALSE = logical(0);

if nargin<2 CONTENTS = FALSE; end;

TESTREGEX=logical(0);

% DEFINE REGEX FOR FUNCTION NAMES
% This thing's pretty ugly but it seems to work.

% GET INPUT
if ~TESTREGEX
    if CONTENTS
        slurp = FIDSlurp;
    else
        slurp = fscanf(FIDSlurp,'%c');
    end;
else
    slurp = { ...
    sprintf(['function [a,b,c,d] = Abc321_a(a,b)\n' ...
    'function [out1, out2] = funname1(in1, in2)\n']), ...
    sprintf(['function funname2(in1, in2)\n' ...
    'function funname3\n']), ...
    sprintf(['function x12 = funname4(in1, in2)\n' ...
    'function A_8_= funname5\n' ...
    'function x12 = funname4(in1, in2)\n' ...
    'function A8_=funname5\n'])};
end;               

% A subfunction is created by defining a new
% function with the function keyword after the body
% of the preceding function or subfunction.
rs.key = 'function\s+(?# line must start with "function")';

% Function, subfunction and variable names must begin
% with a letter, which may be followed by any
% combination of letters, digits, and underscores.
rs.var = '\s*[a-zA-Z]\w*\s*(?# valid variable)' ;

% Output parameters may be of the form "a = function"
% << PLEASE EXCUSE THIS KLUDGE >> Time is money, put a bag
% on the side of it and get it out the door!
slurp = regexprep(slurp,'(function\s+)([^=\[\]]+)=','$1[$2]=');
% or "[a,b,c] = function"
rs.listout = [ '\[' rs.var '(,' rs.var ')*\]' ...
           '(?# a bracketed series of variable names)'];
% combining "a=" with "[a,b,c]="
rs.out = [rs.listout '\s*=(?# fun outputs)'];

% put it all together
rs.prefun = ['(' rs.key '(' rs.out ')+|(' rs.key '))(?# output assignment)'];
rs.funname = ['(?<=' rs.prefun ')' rs.var];

foo = regexp(slurp,rs.funname,'match');

% << PLEASE EXCUSE THIS KLUDGE >> Time is money, put a bag
% on the side of it and get it out the door!
if ~iscellstr(foo)
    foofoo={};
    for k=1:length(slurp)
        foofoo = {foofoo{:} foo{k}{:}};
    end;
foo = foofoo;
end;
foo = regexprep(foo,'(\s+)','');
