%% Make a dictionary of shape descriptors

%%
clc;
clear all;
close all;
addpath('../../');

% Path to image file of letters
imagepath = '~/Data/images/Letters';

nbPoints = 2^10; % Number of boundary points for discretization
delta = 1; % standard size
ord = 5; % order of CGPT dictionary

mode = 1; % 0 for small dico, 1 for dico of letters

%%
% All shapes have the same conductivity and permittivity values
cnd = 5;
pmtt = 2;
lambda = asymp.CGPT.lambda(cnd, pmtt); % contrast constant

%% Make dictionary

disp('Construction of the dictionary...');

B={};

if mode == 0
    B{1} = shape.Ellipse(delta/2,delta/2,nbPoints); % disk
    B{2} = shape.Ellipse(delta*1,delta/2,nbPoints); % ellipse
    B{3} = shape.Flower(delta/2,delta/2,nbPoints); % flower
    B{4} = shape.Triangle(delta, pi/3, nbPoints); % triangle
    B{5} = shape.Rectangle(delta, delta, nbPoints); % square
    B{6} = shape.Rectangle(delta/2, delta, nbPoints); % rectangle
    B{7} = shape.Imgshape([imagepath,'/A.png'], nbPoints); % A
    B{8} = shape.Imgshape([imagepath,'/E.png'], nbPoints); % E
else
    % Or a dictionary of 26 letters
    letters='ABCDEFGHIJKLMNOPQRSTUVWXYZ';

    for n=1:26
        B{n} = shape.Imgshape([imagepath,'/', letters(n),'.png'], nbPoints);
    end
end

%%
% Names of dictionary elements
names = {};
for n=1:lendico
    names{n} = B{n}.name_str;
end

%%
% Calculate the CGPT dictionary
CGPT={}; I1={}; I2={};

for n = 1:length(B)
    CGPT{n} = asymp.CGPT.theoretical_CGPT(B{n}, lambda, ord);

    [I1{n}, I2{n}, ~] = dico.CGPT.ShapeDescriptor_CGPT(CGPT{n});    
end

%% 
% Draw the dictionary shapes and save the dictionary

Dico.ord = ord;
Dico.B = B;
Dico.cnd = cnd;
Dico.pmtt = pmtt;
Dico.lambda = lambda;
Dico.CGPT = CGPT;
Dico.I1 = I1;
Dico.I2 = I2;
Dico.names = names;

if mode == 0
    % for small dico
    for n=1:length(B)
        fig=figure; plot(B{n}, 'LineWidth', 2); axis image;
        saveas(fig, ['~/Data/outputs/figures/small_dico/',B{n}.name_str,'.eps'], 'psc2');    
    end
    save('~/Data/dico/CGPT/smalldico_CGPT.mat','Dico','-v7.3');
else
    % for letters
    for n=1:length(B)
        fig=figure; plot(B{n}, 'LineWidth', 2); axis image;
        saveas(fig, ['~/Data/outputs/figures/letter_dico/',B{n}.name_str,'.eps'], 'psc2');
    end
    save('~/Data/dico/CGPT/letterdico_CGPT.mat','Dico','-v7.3');
end
