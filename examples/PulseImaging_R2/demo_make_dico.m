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
% B{2} = shape.Flower(delta/2,delta/2,nbPoints); % flower
% B{3} = shape.Triangle(delta, pi/3, nbPoints); % triangle
% B{4} = shape.Rectangle(delta/3, delta, nbPoints); % rectangle
% 
% % All shapes have the same conductivity and permittivity values
% % Sea water: 5 Siemens/meter, Fish: 10^-2 S/m
% cnd = [10*ones(1,4)];
% % cnd = [10*ones(1,8), 5];
% pmtt = [ones(1,4)]; % This value is not fixed absolutely. Interaction with the frequency.

% % Full dictionary
B{1} = shape.Ellipse(delta/2,delta/2,nbPoints); % disk
B{2} = shape.Ellipse(delta*1,delta/2,nbPoints); % ellipse
B{3} = shape.Flower(delta/2,delta/2,nbPoints); % flower
B{4} = shape.Triangle(delta, pi/3, nbPoints); % triangle
B{5} = shape.Rectangle(delta, delta, nbPoints); % square
B{6} = shape.Rectangle(delta/3, delta, nbPoints); % rectangle
B{7} = shape.Rectangle(delta/2, delta, nbPoints); % rectangle 2
B{8} = shape.Imgshape([imagepath,'/A.png'], nbPoints); % A

% B{9} = shape.Imgshape([imagepath,'/B.png'], nbPoints); % B
B{9} = shape.Imgshape([imagepath,'/E.png'], nbPoints); % E
B{10} = shape.Imgshape([imagepath,'/L.png'], nbPoints); % L
B{11} = shape.Ellipse(delta*1,delta/2,nbPoints); % ellipse 2 with different cnd and pmtt values

% All shapes have the same conductivity and permittivity values
% Sea water: 5 Siemens/meter, Fish: 10^-2 S/m
% cnd = [10^-2*ones(1,9), 10^-1];
cnd = [10*ones(1,length(B)-1), 5];
pmtt = [ones(1,length(B)-1), 2]; % This value is not fixed absolutely. Interaction with the frequency.

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

% Scl = 1.5.^(-6:-1); % Best scales found by data analysis for cnd=[10^-2...10^-2,10^-1], pmtt=[1..1,2]
% Scl = 2.^(-3:2); % Best scales found by data analysis for cnd=[2...2,5], pmtt=[1..1,2]. Larger than 2, the ellipse is shrinked to a circle, smaller than -3, rotational symmetric objects are all similar to each other.
% Scl = 2.^(-1:2); % Best scales found by data analysis for cnd=[10...10,5], pmtt=[1..1,2]
Scl = 2.^(-2:3);
nbScl = length(Scl); % number of scales

Ntime = 2^10; % time interval length 2^10
Tmax0 = 5;
Tmax = zeros(1, nbScl);
dt0 = zeros(1, nbScl);

Fmax = zeros(1, nbScl);
df = zeros(1, nbScl);

waveform = zeros(nbScl, Ntime);
freqform = zeros(nbScl, Ntime);

CGPTt0 = cell(length(B), nbScl);
SDt = cell(length(B),1);  % Shape descriptor

extrema = {};

for s = 1:nbScl
    % pulse waveform at the scale s    
    [waveform(s,:), dt(s), Tmax(s), freqform(s, :), df(s), Fmax(s), extrema{s}] = tools.make_pulse(Tmax0, Ntime, Scl(s));
end

% Compute CGPT
for m=1:length(B) % iteration on the shape
    
    fprintf('Proceeding the shape %s...\n', B{m}.name_str);
    
    % Compute time-dependent CGPT
    
    tic
    for s = 1:nbScl
        fprintf('...Proceeding the scale %f...\n', Scl(s));
        [CGPTt0{m,s}, dt0(s), ~] = asymp.CGPT.theoretical_CGPT_time(B{m}, cnd(m), pmtt(m), ord, ...
                                                          freqform(s,:), df(s)); % dt depends only on df
        
        % % Compute CGPT by reconstruction, faster but less accurate
        % [CGPTt0{m,s}, dt0(s), ~] = asymp.CGPT.theoretical_CGPT_time_recon(B{m}, cnd(m),
        % pmtt(m), ord, Tmax0, Ntime, Scl(s));
        
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
% Dico.comments = 'The value of conductivities is based on real values of sea water and fish. The value of scaling is based on numerical tuning: it is the range on which the similar shapes are most distinct.';

fname = ['/Volumes/ExFAT200G/Data/dico/Pulse/smalldico',num2str(length(Dico.B)),'_', num2str(nbScl),'scl.mat'];
save(fname,'Dico','-v7.3');
fprintf('Data saved in %s\n', fname);

%% manipulations
% pathname = '/Volumes/Macbook/Data/';
% load([pathname,'/dico/Pulse/smalldico9_21scl.mat']);
% 
% Dico0 = Dico;
% sidx = 10:13;
% nbScl = length(sidx);
% 
% Dico.Scl = Dico0.Scl(sidx);
% Dico.Tmax = Dico0.Tmax(sidx);
% Dico.dt0 = Dico0.dt0(sidx);
% Dico.Fmax = Dico0.Fmax(sidx);
% Dico.df = Dico0.df(sidx);
% Dico.waveform = Dico0.waveform(sidx,:);
% Dico.extrema = Dico0.extrema(sidx);
% Dico.freqform = Dico0.freqform(sidx, :);
% Dico.CGPTt0 = Dico0.CGPTt0(:,sidx);
% 
% fname = ['/Volumes/Macbook/Data/dico/Pulse/smalldico',num2str(length(Dico.B)),'_', num2str(nbScl),'scl.mat'];
% save(fname,'Dico','-v7.3');
% fprintf('Data saved in %s\n', fname);
% 
