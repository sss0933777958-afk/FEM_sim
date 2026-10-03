%% xa_cross.m -- is R=450 below R=350 because of the DATA or because of the KERNEL?
%  2x2: each of the two ranges fitted with each of the two winning kernels.  If the
%  peak tracks the kernel, the difference is regularisation; if it tracks the range,
%  it is the data.  [ADDED 2026-09-16, diagnostic]
clear; clc;
here=fileparts(mfilename('fullpath')); FMX=fileparts(fileparts(here));
MAIN=fileparts(fileparts(fileparts(FMX))); FLX=fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
Pq = [linspace(-150,150,601).' zeros(601,2)];
I  = zeros(6,1); I(1) = 1;
RR = [350 450];  KK = [500 1e-3; 610 1e-2];
fprintf(['%6s %8s %9s %10s %12s' newline],'R','rho','lam','nodes','peak(+150)');
for a = 1:2
  for b = 1:2
    [~,gq,iq] = rbf_field(RR(a), struct('rho',KK(b,1),'lam',KK(b,2)));
    G = gq(Pq, I, 'du');
    fprintf(['%6d %8g %9.0e %10d %12.4f' newline], RR(a), KK(b,1), KK(b,2), size(iq.P,1), G(end,1));
    clear gq iq
  end
end
