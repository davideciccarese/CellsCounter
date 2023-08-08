%% Davide Ciccarese

% Date of creatinon: 12/07/2023
% Last Modification: 08/08/2023

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
cd 'C:\Users\Van der Meer group\Desktop\'


%% FIRST STEP parameters

%Parametrization following Shaerly protocol
lambda = 2; % IMPORTANT LAMBDA IS THE AVARAGE number of cells/droplet; Schaerli used 0.3 
r = 20*10^-6; %25*10^-6 Um, radius of Beads

%FACS cells in 1 ml
ODe1 = 3*10^8; %NatCom
ODe2 = 2*10^9; %Fluorecence species


%%

vecT = [ODe1,ODe2];
% ODei  = (vecT.*([vecT==(min(vecT))]));
% ODe  = ODei(ODei>0);

varMix = {'Dilute NatCom ', 'Dilute Fluo '};
% disp(varMix(([vecT==(min(vecT))])))

Sel1 = (vecT.*([vecT==(min(vecT))]));
Sel1 = Sel1(Sel1>0);
Sel2 = (vecT.*([vecT==(max(vecT))]));
Sel2 = Sel2(Sel2>0);

Dilution = round((Sel1/Sel2)*2000);
Dilution2 = 2000-Dilution;

% Poisson distribution

lambda % Parameter for the Poisson distribution

%----EXAMPLE----
% ? = 0.3 ? P(0, 0.3) = 0.74; P(1, 0.3) = 0.22; P(2, 0.3) = 0.033; P(3, 0.3) = 0.003.
% This means, if on average we have 0.3 cells/droplet we will have 74%
% empty droplets, 22% droplets with 1 cell, 3.3% droplets with 2 cells and
% 0.3% droplets with 3 cells.

k = 0:10;     % Values at which to evaluate the distribution number of cells in a droplet

% Calculate the probability mass function (PMF) for the Poisson distribution
pmf = poisspdf(k, lambda);

% Display the results
% disp('Poisson PMF:');
Percent = round(pmf*100,2); %in percentage
numbCells = k; % Values at which to evaluate the distribution number of cells in a droplet


% Visualize the distribution
bar(k,Percent,1);
saveas(gcf,'Distribution.png')
 
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
writetable(T, 'Cells_number_distribution.txt','Delimiter', '\t');

X = cell2mat([((varMix(([vecT==(max(vecT))])))),'',num2str(Dilution),' ul in ',num2str(Dilution2),' ul of PBS'])

% varMs = {'SymCom', 'Fluo'};
% %varM = [(varMs(([vecT==(min(vecT))]))), 'and'  (varMs(([vecT==(max(vecT))]))), '1:1']
varM = ['WARNING!!....Then mix NatCom & Fluo 1:1, e.g. 500 ul and 500 ul and measure again FACS change mOD RUN second step ']



fNlMessage = {char(X); varM}

fileID = fopen('Cells_number_distribution.txt', 'a');
for i = 1:numel(fNlMessage)
    fprintf(fileID, '%s\n', fNlMessage{i});
end
fclose(fileID);

% "These calculations will give you a first indication for your cell
% density. If the flow cytometry analysis indicates that too
% many or too few beads contain cells, the cell density should be
% accordingly adjusted in the next experiment."



%% Second steps

%================================
% FACS after mixing 1:1
lambda = 2; % IMPORTANT LAMBDA IS THE AVARAGE number of cells/droplet; Schaerli used 0.3 
r = 20*10^-6; %25*10^-6 Um, radius of Beads
agarosePerc = 2;
mOD = 2*10^8;
FinalVol = 3000;
%%



% Estimated OD to cell per droplets

V = 4/3*(pi*r^3); %volume droplet
L = V*10^3; %Liter volume droplets

dL = 1/L; % n. Droplet per liter
dLml = dL/10^3; %n. Droplet per milliliter

TotNCells = lambda*dLml; %how many cell per Droplet

cOD = TotNCells/(mOD); %Dilution factor


%---How many time you need to dilue your cells!!---
% (mOD/cOD)
% foldDilution = round(,3);
foldDilution = cOD;
Dilution3 = round(FinalVol*foldDilution);


DivFactAg = 1/agarosePerc; %dilution factor of 2% agarose vs 1% agarose
AgC = round((DivFactAg*FinalVol));

pbsAgar = FinalVol-(AgC+Dilution3);

% Dilution4 = 1000-Dilution3;

finalMessage = ['Dilute your 1:1 mix ', num2str(Dilution3), ' ul in ', num2str(AgC), 'uL of ',num2str(agarosePerc), '% Agarose', ' then add ',num2str(pbsAgar),'uL of PBS', ' For a final volume of ',num2str(FinalVol)]


