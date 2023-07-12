%% Davide Ciccarese

% Date of creatinon: 12/07/2023
% Last Modification: 12/07/2023

% The code calculates the required dilution factor for a given cell
% density, based on parameters such as the average number of cells per
% droplet and the measured optical density. It also calculates the
% probability of having a certain number of cells in a droplet using a
% Poisson distribution. The code, is based on the assumption of a typical
% OD of E.coli. This can be changed using FACS
% 
% ---Reference---
% The code is basd on this paper
% Duarte JM, Barbier I, Schaerli Y. (2017) Bacterial Microcolonies in Gel Beads
% for High-Throughput Screening of Libraries in Synthetic Biology. ACS
% Synth Biol; 6: 1988?1995.

%% Parameters

clear all
close all


%Parameters from Shaerly protocol

ODe = 5*10^8;%typical E.coli density OD 600= 1  is 5 x 10^8 cells/ml.
lambda = 0.3 % Shaerly used 0.3 IMPORTANT lambda avarage n of cells/droplet
r = 25*10^-6; %?m, radius of droplets
mOD = 0.2; %measured OD;

%% Poisson distribution

lambda % Parameter for the Poisson distribution

%----EXAMPLE----
% ? = 0.3 ? P(0, 0.3) = 0.74; P(1, 0.3) = 0.22; P(2, 0.3) = 0.033; P(3, 0.3) = 0.003.
% This means, if on average we have 0.3 cells/droplet we will have 74%
% empty droplets, 22% droplets with 1 cell, 3.3% droplets with 2 cells and
% 0.3% droplets with 3 cells.

k = 0:10;     % Values at which to evaluate the distribution number of cells in a droplet

% Calculate the probability mass function (PMF) for the Poisson distribution
pmf = poisspdf(k, lambda);

% Visualize the distribution
bar(k,Percent,1);
saveas(gcf,'Distribution.png')

% Display the results
% disp('Poisson PMF:');
Percent = (pmf*100); %in percentage
numbCells = 0:length(pmf)-1;

tbl = [Percent numbCells];
% % Specify variable names
% VariableNames = {'Percentage', 'Number_of_cells_per_beads'};
% 
% % Create the table
% tbl = table(Percent, numbCells);
% % Create the table
% tbl = table(Percent, numbCells, 'VariableNames', VariableNames);

% Create a table with the data and variable names
T = table(Percent', numbCells', 'VariableNames',  {'Percentage', 'Number_of_cells_per_beads'})
% Write data to text file
writetable(T, 'Cells_numbber_distribution.txt','Delimiter', '\t');

%% Estimated OD to cell per droplets

V = 4/3*(pi*r^3); %volume droplet
L = V*10^3; %Liter volume droplets

dL = 1/L; % n. Droplet per liter
dLml = dL/10^3; %n. Droplet per milliliter

TotNCells = lambda*dLml; %how many cell per Droplet

cOD = TotNCells/(ODe);

%---How many time you need to dilue your cells!!---
foldDilution = mOD/cOD

% "These calculations will give you a first indication for your cell
% density. If the flow cytometry analysis indicates that too
% many or too few beads contain cells, the cell density should be
% accordingly adjusted in the next experiment."

