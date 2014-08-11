%% Make a dictionary of shapes and time dependent CGPTs with the class |PulseImaging|

% We need to at least fix one group of values (e.g., permittivity and conductivity) 
% so that the best window of observation in the frequency domain (or equivalently 
% the scaling of the waveform) can be determined numerically. One criteria
% for this is to use similar shapes (eg, ellipse and rectangle) and see on which
% band of frequency these shapes are most easily distinct. From a
% biological point of view, one could think that the electric fish tunes its
% frequency to focus on certain species of preys.

%%
clear all;
close all;
addpath('~/SIES');

% Path to image file of letters
imagepath = '~/Data/images/Letters';

%% Make dictionary

disp('Construction of the dictionary...');

nbPoints = 2^10; % Number of boundary points for discretization 2^10
delta = 1; % standard size

% B{1} = shape.Ellipse(delta/2,delta/2,nbPoints); % disk
% B{2} = shape.Ellipse(delta*1,delta/2,nbPoints); % ellipse
% B{3} = shape.Rectangle(delta/2, delta, nbPoints); % rectangle

% B{1} = shape.Ellipse(delta/2,delta/2,nbPoints); % disk
% B{2} = shape.Flower(delta/2,delta/2,nbPoints); % flower
% B{3} = shape.Triangle(delta, pi/3, nbPoints); % triangle
% B{4} = shape.Rectangle(delta, delta, nbPoints); % square

B{1} = shape.Ellipse(delta/2,delta/2,nbPoints); % disk
B{2} = shape.Ellipse(delta*1,delta/2,nbPoints); % ellipse
B{3} = shape.Flower(delta/2,delta/2,nbPoints); % flower
B{4} = shape.Triangle(delta, pi/3, nbPoints); % triangle
B{5} = shape.Rectangle(delta, delta, nbPoints); % square
B{6} = shape.Rectangle(delta/2, delta, nbPoints); % rectangle
B{7} = shape.Imgshape([imagepath,'/A.png'], nbPoints); % A
B{8} = shape.Imgshape([imagepath,'/E.png'], nbPoints); % E
B{9} = shape.Ellipse(delta*1,delta/2,nbPoints); % ellipse 2 with different cnd and pmtt values

% All shapes have the same conductivity and permittivity values
% Sea water: 5 Siemens/meter, Fish: 10^-2 S/m
cnd = [10^-2*ones(1,8), 10^-1];
pmtt = [ones(1,8), 2]; % This value is not fixed absolutely. Interaction with the frequency.


%%
% Names of dictionary elements
names = {};
for n=1:length(B)
    names{n} = B{n}.name_str;
end

%% Compute multiscale time-dependent CGPT

disp('Computation of theoretical time dependent CGPTs...');
%% 
% Parameters
ord = 1; % order of CGPT dictionary

Scl = 1.5.^(-6:-1);
scl = length(Scl); % number of scales

Ntime = 2^10; % time interval length 2^10
Tmax0 = 5;
Tmax = zeros(1, scl);
dt0 = zeros(1, scl);

Fmax = zeros(1, scl);
df = zeros(1, scl);

waveform = zeros(scl, Ntime);
freqform = zeros(scl, Ntime);

CGPTt0 = cell(length(B), scl);
SDt = cell(length(B),1);  % Shape descriptor

extrema = {};

for s = 1:scl
    % pulse waveform at the scale s    
    [waveform(s,:), dt(s), Tmax(s), freqform(s, :), df(s), Fmax(s), extrema{s}] = tools.make_pulse(Tmax0, Ntime, Scl(s));
end

% Compute CGPT
for m=1:length(B) % iteration on the shape
    
    fprintf('Proceeding the shape %s...\n', B{m}.name_str);
    
    % Compute time-dependent CGPT
    
    tic
    for s = 1:scl
        [CGPTt0{m,s}, dt0(s), ~] = asymp.CGPT.theoretical_CGPT_time(B{m}, cnd(m), pmtt(m), ord, freqform(s,:), df(s)); % dt depends only on df
        
        % [CGPTt{m,s}, dt(s)] = asymp.CGPT.CGPT_time_truncation(CGPTt0, dt0(s), Tmax, Ntime);
    end
    toc
end

Dico = [];

Dico.ord = ord;
Dico.Scl = Scl;
Dico.cnd = cnd;
Dico.pmtt = pmtt;
Dico.B = B;
Dico.names = names;
Dico.Tmax0 = Tmax0;
Dico.Tmax = Tmax;
Dico.dt0 = dt0;
Dico.Fmax = Fmax;
Dico.df = df;
Dico.waveform = waveform;
Dico.extrema = extrema;
Dico.freqform = freqform;
Dico.Ntime = Ntime;
Dico.CGPTt0 = CGPTt0; % The time dependent CGPT in the highest resolution
Dico.comments = 'The value of conductivities is based on real values of sea water and fish. The value of scaling 1.5^(-6:-1) is based on numerical tuning: it is the range on which the ellipse and the rectangle are most distinct.';

fname = ['~/Data/dico/Pulse/smalldico',num2str(length(Dico.B)),'_', num2str(scl),'scl.mat'];
save(fname,'Dico','-v7.3');
fprintf('Data saved in %s\n', fname);

