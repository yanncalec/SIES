function [err, idx, scl] = SCT_matching(S, sfrange, D, frange)
% Dictionary matching algorithm by shape descriptor of SCT
% INPUTS:
% S: shape-descriptor extracted from multi-frequency data, a 3D array.
% sfrange: scanning frequency range of S (used for data acquisition) [sf1, sf2]
% D: the dictionary of pre-computed frequency dependent shape descriptor (inv. to translation and
% rotation), a cell. D{n} is the shape descriptor of the n-th element, a 3d array, and the third
% dimension corresponds to the frenquency.
% frange: frequency range of D [f1, f2]
% OUTPUTS:
% err: euclidean distance between S and each element of D
% idx: sorted dictionary elements in decreasing order of similarity

err = zeros(1, length(D));
scl = zeros(1, length(D));

for n=1:length(D)
    mD = squeeze(mean(mean(D{n},1),2));
    mS = squeeze(mean(mean(S,1),2));

    [serr, vscl, idx] = dico.Helmholtz.scaling_lookup_table(mD, frange, mS, sfrange);
    err(n) = serr(idx);
    scl(n) = vscl(idx);
end

[~,idx]=sort(err);

