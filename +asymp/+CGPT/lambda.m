function val = lambda(cnd, pmtt, freq)
% Compute the contrast lambda from the conductivity, the permittivity and the frequency.

if nargin < 3
	freq = 0;
end

if nargin < 2 || isempty(pmtt)
	pmtt=zeros(size(cnd));
end

if ~isscalar(freq) || ~isreal(freq) || freq<0
	error('Frequency must be a positive scalar.');
end

if iscell(cnd)
	cnd=cell2mat(cnd);
end

for n=1:length(cnd)
	if cnd(n) == 1 || cnd(n)<0
		error('Invalid value of conductivity.');
	end
end

%%%%%%%%%%%%%%%%%%%%%%%% Important Remark %%%%%%%%%%%%%%%%%%%%%%%
% The facter 2*pi comes from the Fourier transform's convention:
%   f^(w) = \int f(x) exp(-2*pi*1i*x*w) dx    (1)
% If the convention
%   f^(w) = \int f(x) exp(-1i*x*w) dx         (2)
% is used, then there will be no 2*pi factor.
% This affects notably the class PulseImaging_R2 (though the shape
% identification with multifreqeuncy data could be affected too). In
% published papers the convention (2) is used, however in the library
% it is the convention (1) which is implemented (to be compatible with the
% fft of matlab)

toto = cnd + 2*pi*1i*pmtt*freq; % convention (1)
% toto = cnd + 1i*pmtt*freq; % convention (2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

val = (toto+1)./(toto-1)/2;