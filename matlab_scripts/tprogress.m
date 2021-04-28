
function tprogress(p,varargin)
%TPROGRESS shows percent completed on screen
% TPROGRESS(p) is a text version of "waitbar" and similar functions. 
% It shows the percent completed of a process directly on 
% the screen as a %. 'p' must be a proportion (0 to 1). A text string 
% can be supplied as optional argument. It can handle multiple counters,
% just make 'p' a vector and the optional text string (if used) a cell
% array. Make sure 'p' and the cell array have the same number of
% elements.
% All 'p' should start at 0 and end at 1 for it to work correctly.
% This a minimalistic function with simplicity in mind.
% Not many options, no argument check. No initialization necessary.
% 
% Example:
% M= 5; N= 25;
% for m= 0:M
% for n= 0:N
%		pause(0.02) % substitute here your code
%		tprogress([n/N,m/M],{'First ';'Secnd '})
% end
% end
% 
%  Author: Francisco de Castro. Aug. 2012
%	Inspired by: tprogressbar, Paul Proteus
% Default prompts
text = repmat({'Prog: '},numel(p),1);
% User-supplied prompt
if ~isempty(varargin) text= varargin{1}; end
if ~iscell(text) text= {text}; end
	
% Put p in percent. Problems if p > 1. Un-comment next line if you want a
% warning about this.
% if p > 1 warning('Argument is > 1'); end
p = 100* abs(p); 
% Initialise string
strpcent= []; 
% Fill-in the string
for j= 1:numel(p)
	strpcent= [strpcent text{j} num2str(p(j),'%05.1f') '%% '];
end
backsp = repmat('\b',1,length(strpcent)-numel(p));
whites = repmat(' ', 1,length(strpcent)-numel(p));
if sum(p) == 0 fprintf(1,whites); end
fprintf(1,[backsp strpcent])
if all(p == 100) fprintf(1,'Done\n'), end
