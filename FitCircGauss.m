function [optparam, fo, circGauss] = FitCircGauss(x, y)
% [optparam, fo, sinfunc] = FitCircGauss(x, y)
%
% 
% Fitting Data to following function with "lsqnonlin"
% circGauss = @(v,x) v(1).*exp(v(4)*(cos(2.*(x-v(3)))-1)) + v(2);   
%
%   x: xaxis (radian)
%   y: actuial value
%   
%   optparam: optimal parameter(Amp, offset, phase)
%   fo : norm^2 
%   circGauss: function handle of circGauss(param, x) 
%
% example codes
%
% yfunc  = @(v,x) v(1).*exp(v(4)*(cos(2.*(x-v(3)))-1)) + v(2);   
% yparam = [10,1,pi/4,pi/3]; 
% x = pi*(0:15:180)/180;
% y = randn(size(x)) + yfunc(yparam, x);
% 
% [optparam, ~, circGauss]=FitCircGauss(x,y);
% xopt = (0:1:180)*pi/180;
% figure; hold on
% plot(x,y,'d')
% plot(xopt, circGauss(optparam, xopt), 'k-')
%
% 2016-01-21 Ryosuke F Takeuchi

circGauss = @(v,x) v(1).*exp(v(4)*(cos(2.*(x-v(3)))-1)) + v(2);   

x = x(:);
y = y(:);
[~,init_phase] = max(y);
vinit(1) = max(y)-min(y); % Amplitude
vinit(2) = min(y);        % offset
vinit(3) = x(init_phase);   % x-phase
vinit(4) = 0.5*pi;        % sigma of gaussian

options = optimoptions('lsqnonlin', 'MaxIter', 10000,...
	'Algorithm', 'trust-region-reflective', 'Tolx', 10^-15, 'TolFun', 10^-15,...
	'ScaleProblem', 'none', 'MaxFunEvals', length(x)*200, 'Display', 'off');

ub = [abs(vinit(1)*2), max(y), 2*vinit(1) 180];
lb = [0, -abs(max(y)), 0 0];
	
initphz = (0:30:360)*pi/180;
efunc = @(v) abs((circGauss(v,x) - y));

for i = 1:length(initphz)
	vinit_tmp = vinit;
	vinit_tmp(3) = initphz(i);
	[params{i}, f(i), flag{i}, lamb{i}] = ...
      lsqnonlin(efunc, vinit_tmp, lb, ub,options);	
end

[~, o] = min(f);
optparam = params{o};
fo = f(o);

