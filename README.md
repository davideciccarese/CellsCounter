# CellsCounter for droplet preparation
'CellsCalculator' code calculates the required dilution factor for a given cell density, based on parameters such as the average number of cells per droplet and the measured optical density. It then calculates the probability of having a certain number of cells in a droplet based on the Poisson distribution. The code, is based on the assumption of a typical number of cells of E.coli per ml. But can be changed using FACS.

'CellsCalculator_2Species_FACS' code calculates the necessary dilution factor for a given cell density, which is measured using FACS, for a mixture of two species. It then calculates the probability of having a certain number of cells in a droplet based on the Poisson distribution.

'CellsCalculator_2Species' code calculates the necessary dilution factor for a given cell density, which is measured using optical density (OD 600), for a mixture of two species. It then calculates the probability of having a certain number of cells in a droplet based on the Poisson distribution. 

---Reference---
These codes are basd on this paper
Duarte JM, Barbier I, Schaerli Y. (2017) Bacterial Microcolonies in Gel Beads
for High-Throughput Screening of Libraries in Synthetic Biology. ACS
Synth Biol; 6: 1988–1995.

