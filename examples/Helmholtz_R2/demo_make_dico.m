%% Construction of SCT dictionary
% Make a small dictionary of shape descriptors using multiple frequencies.

%%
clc;
clear all;
close all;
addpath('../../');

% Path to image file of letters
imagepath = '~/Data/images/Letters';

nbPoints = 2^10; % Number of boundary points for discretization
delta = 1; % standard size

%%
% All shapes have the same permeability and permittivity values
pmtt = 3; pmeb = 3;

%%
% Background values
pmtt_bg = 1; pmeb_bg = 1;

%% Make dictionary

disp('Construction of the dictionary...');

B={};

B{1} = shape.Ellipse(delta/2,delta/2,nbPoints); % disk
B{2} = shape.Ellipse(delta*1,delta/2,nbPoints); % ellipse
B{3} = shape.Flower(delta/2,delta/2,nbPoints); % flower
B{4} = shape.Triangle(delta, pi/3, nbPoints); % triangle
B{5} = shape.Rectangle(delta, delta, nbPoints); % square
B{6} = shape.Rectangle(delta/2, delta, nbPoints); % rectangle
B{7} = shape.Imgshape([imagepath,'/A.png'], nbPoints); % A
B{8} = shape.Imgshape([imagepath,'/E.png'], nbPoints); % E

%%
% Names of dictionary elements
names = {};
for n=1:length(B)
    names{n} = B{n}.name_str;
end

%% Compute the SCT dictionary

%%
% Apriori range of scaling
sclrange = [1/2, 2];

%% 
% Range of scanning frequencies
sfrange = [1, 2]*pi;

%% 
% Frequency range of the dictionary is deduced from the apriori range of
% the scaling and the range of scanning requency
frange = sclrange .* sfrange;
Nf = floor((frange(2)-frange(1))/0.05);  % number of dictionary frequency
freq = linspace(frange(1), frange(2), Nf);

%%
% Order of the dictionary
ord = 30;

%% 
% dimension of the shape descriptor
Nv = 512; 

%%
% Compute the SCT and shape descriptor
for n=1:length(B)
    fprintf('Processing the shape %d...\n', n);

    % For each frequency
    SCT = zeros(2*ord+1, 2*ord+1, Nf);
    for f=1:Nf
        SCT(:,:,f) = asymp.SCT.theoretical_SCT(B{n}, pmeb, pmtt, pmeb_bg, pmtt_bg, ord, freq(f));
    end
    
    [S, G] = dico.SCT.ShapeDescriptor_SCT(SCT, Nv);

    Dico.SCT{n} = SCT; % SCT
    Dico.SD_S{n} = S; % Shape descriptor S
    Dico.SD_G{n} = G; % Shape descriptor G    
end

clear SCT S G;

Dico.sfrange = sfrange;
Dico.sclrange = sclrange;
Dico.frange = frange;
Dico.freq = freq;

Dico.ord = ord;
Dico,Nv = Nv;
Dico.B = B;
Dico.pmtt = pmtt;
Dico.pmeb = pmeb;
Dico.pmtt_bg = pmtt_bg;
Dico.pmeb_bg = pmeb_bg;
Dico.names = names;

save('~/Data/dico/SCT/smalldico_SCT.mat','Dico','-v7.3');

