function val = wavenb(freq, pmeb, pmtt)
% Wavenumber
% pmeb, pmtt: permeability and permittivity
    val = freq*sqrt(pmtt*pmeb);
end