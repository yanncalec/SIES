clc;
clear all;
close all;
addpath('~/OOP/');
imagepath = '~/Data/images/Letters/';

%% Construction of the dictionary
disp('Construction of the dictionary...');

nbPoints = 2^10;
delta = 1;

D={};

D{1} = shape.Ellipse(delta/2,delta/2,[0,0]',0,nbPoints); % disk
D{2} = shape.Triangle(delta, pi/3, nbPoints);
D{3} = shape.Rectangle(delta, delta, nbPoints); % square
D{4} = shape.Ellipse(delta*1,delta/2,[0,0]',0,nbPoints); % ellipse
D{5} = shape.Flower(delta/2,delta/2,[0,0]',0,nbPoints,5,0.4,0);
D{6} = shape.Imgshape([imagepath,'A.png'], delta, delta, nbPoints);
D{7} = shape.Imgshape([imagepath,'E.png'], delta, delta, nbPoints);
D{8} = shape.Rectangle(delta/2, delta, nbPoints);

%% draw the dictionary shapes
% fig = figure;
% subplot(241); plot(D{1}); axis image;
% subplot(242); plot(D{2}); axis image;
% subplot(243); plot(D{6}); axis image;
% subplot(244); plot(D{7}); axis image;
% subplot(245); plot(D{3}); axis image;
% subplot(246); plot(D{8}, '--'); axis image;
% subplot(247); plot(D{4}); axis image; xlim([-delta,delta]/2*1.1); ylim([-0.4,0.6]);
% subplot(248); plot(D{5}); axis image; xlim([-delta,delta]/2*1.1); ylim([-delta,delta]/2*1.1);
% saveas(fig, 'dico.eps', 'psc2');

% draw separately
fig=figure; plot(D{1}); axis image; 
xlim([-delta,delta]*0.55); ylim([[-delta,delta]*0.55]); 
saveas(fig, ['../figures/', D{1}.name_str, '.eps'], 'psc2');

fig=figure; plot(D{2}); axis image; 
xlim([-delta,delta]*0.55); ylim([-0.35, 0.75]); 
saveas(fig, ['../figures/', D{2}.name_str, '.eps'], 'psc2');

fig=figure; plot(D{3}); axis image; 
xlim([-delta,delta]*0.55); ylim([[-delta,delta]*0.55]); 
saveas(fig, ['../figures/', D{3}.name_str, '.eps'], 'psc2');

fig=figure; plot(D{4}); axis image; 
xlim([-delta,delta]*1.05); ylim([[-delta,delta]*1.05]); 
saveas(fig, ['../figures/', D{4}.name_str, '.eps'], 'psc2');

fig=figure; plot(D{5}); axis image; 
xlim([-delta,delta]*0.7); ylim([[-delta,delta]*0.7]); 
saveas(fig, ['../figures/', D{5}.name_str, '.eps'], 'psc2');

fig=figure; plot(D{6}); axis image; 
xlim([-delta,delta]*0.55); ylim([[-delta,delta]*0.55]); 
saveas(fig, ['../figures/', D{6}.name_str, '.eps'], 'psc2');

fig=figure; plot(D{7}); axis image; 
xlim([-delta,delta]*0.55); ylim([[-delta,delta]*0.55]); 
saveas(fig, ['../figures/', D{7}.name_str, '.eps'], 'psc2');

fig=figure; plot(D{8}); axis image; 
xlim([-delta,delta]*0.55); ylim([[-delta,delta]*0.55]); 
saveas(fig, ['../figures/', D{8}.name_str, '.eps'], 'psc2');

% %%  draw separately
% for n=1:length(D)    
%     fig = figure; 
%     plot(D{n}); axis image;
%     xlim([-delta,delta]*0.7); ylim([-delta,delta]*0.7);
%     saveas(fig, ['../figures/', D{n}.name_str, '.eps'], 'psc2');
% end
