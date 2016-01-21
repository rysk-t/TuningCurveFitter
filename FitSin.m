function [optparam, fo, sinfunc] = FitSin(x, y)
% [optparam, fo, sinfunc] = FitSin(x, y)
%
% 
% Fitting Data to following function with "lsqnonlin"
% sinfunc = @(v, x) (v(1)*sin(x-v(2))+v(3));
%
%   x: xaxis (radian)
%   y: actuial value
%   
%   optparam: optimal parameter(Amp, offset, phase)
%   fo : norm^2 
%   sinfunc: function handle of sinfunc(param, x) 
%
% x = pi*(0:45:360)/180;
% y = 2*rand(1)*sin(x-randn(1)*pi+pi/2)+0.2;
% yraw = y+0.6*randn(size(x));
% 
% [optparam, ~, sinfunc]=FitSin(x,y);
% xopt = (0:1:360)*pi/180;
% figure; hold on
% plot(x,yraw,'d')
% plot(xopt, sinfunc(optparam,xopt), 'k-')
%
% 2015-12-09 Ryosuke F Takeuchi

x = x(:);
y = y(:);

sinfunc = @(v, x) (v(1)*sin(x-v(2))+v(3));
vinit(1) = max(y)+mean(y); % Amplitude
vinit(2) = 0;              % offset
vinit(3) = max(y);         % x-phase 

options = optimoptions('lsqnonlin', 'MaxIter', 10000,...
	'Algorithm', 'trust-region-reflective', 'Tolx', 10^-15, 'TolFun', 10^-15,...
	'ScaleProblem', 'none', 'MaxFunEvals', length(x)*200, 'Display', 'off');

ub = [vinit(1)*2, 2*pi, vinit(1)];
lb = [0, 0, min(y)];
	
initphz = (0:30:360)*pi/180;
efunc = @(v) abs((sinfunc(v,x) - y));

for i = 1:length(initphz)
	vinit(2) = initphz(i);
	[params{i}, f(i), flag{i}, lamb{i}] = ...
      lsqnonlin(efunc, vinit, lb, ub,options);	
end

[~, o] = min(f);
optparam = params{o};
fo = f(o);
