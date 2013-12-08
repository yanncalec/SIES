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
    
toto = cnd + 1i*pmtt*freq;
val = (toto+1)./(toto-1)/2;