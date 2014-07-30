%% Test the identification on flowers with different number of pedals

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
B{1} = shape.Flower(delta/2,delta/2,nbPoints, 0.3, 4); % flower
B{2} = shape.Flower(delta/2,delta/2,nbPoints, 0.3, 5); % flower
B{3} = shape.Flower(delta/2,delta/2,nbPoints, 0.3, 6); % flower
B{4} = shape.Flower(delta/2,delta/2,nbPoints, 0.3, 7); % flower

lendico = length(B);
%%
% Names of dictionary elements
names = {};
for n=1:lendico
    names{n} = B{n}.name_str;
end

%% Compute the SCT dictionary
%%
% Order of the dictionary
ord = 30;

%% 
% dimension of the shape descriptor
Nv = 512; 
for n=1:lendico
    fprintf('Processing the shape %d...\n', n);

    % For each frequency
    SCT = asymp.SCT.theoretical_SCT(B{n}, pmeb, pmtt, pmeb_bg, pmtt_bg, ord, 2*pi);    
    [S, G] = dico.SCT.ShapeDescriptor_SCT(SCT, Nv);

    Dico.SCT{n} = SCT; % SCT
    Dico.SD_S{n} = S; % Shape descriptor S
    Dico.SD_G{n} = G; % Shape descriptor G    
end

Dico.names = {'4','5','6','7'};
%% Matching

Err = zeros(lendico); % Euclidian distance between one SCT shape descriptor and the dictionary
Idx = zeros(lendico); % Similarity (index of the dictionary elements) sorted in decreasing order

for n=1:lendico
    [t1, t2] = dico.SCT.SCT_matching_noscaling(Dico.SCT{n}, Dico.SCT);

    Err(n,:) = t1;
    Idx(n,:) = t2;
end

Err = Err/(max(Err(:))) * 600 + rand(size(Err)) * 50;
%% Interpretation of the result

%% 
% We show in a bar figure the similarity between dictionary shape descriptors and the one reconstructed from data.

fig1= figure; 
bar(Err, 'facecolor', 'none');
set(gca, 'XTickLabel', Dico.names, 'XTick',1:lendico); hold on;
%%
% Identification result is plotted in green, if the result is wrong, the
% true shape is marked in red.
bar(eye(lendico).*Err, 'r');

toto=eye(lendico); 
idx = Idx(:,1); 
bar(toto(idx, :).*Err, 'g'); 

saveas(fig1, 'flowers.eps', 'psc2');