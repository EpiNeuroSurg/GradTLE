% load data
Demo_HCs = readtable('/host/verges/tank/data/mo/BIDS_CERMEP-MRXFDG/BIDS_CERMEP-MRXFDG_demo.csv');
Demo_disPTs = readtable('/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/dis_BTTH1_BFTH.csv');
Demo_valPTs = readtable('/host/verges/tank/data/mo/BIDS_20260316_TLE_BTTH2_BSBH/val_BTTH2_BSBH.csv');
Demo_senPTs = readtable('/host/verges/tank/data/mo/BIDS_20260531_TLE_BTTH3/ses_BTTH3.csv');

Age_HCs     = Demo_HCs.Age; 
Age_disPTs  = Demo_disPTs.AgeSurg;
Age_valPTs  = Demo_valPTs.AgeSurg;
Age_sesPTs  = Demo_senPTs.AgeSurg;
Age_all2    = [Age_HCs; Age_disPTs];
Age_all4    = [Age_HCs; Age_disPTs; Age_valPTs; Age_sesPTs];

Sex_HCs     = Demo_HCs.Sex; 
Sex_disPTs  = Demo_disPTs.Sex;
Sex_valPTs  = Demo_valPTs.Sex; 
Sex_sesPTs  = Demo_senPTs.Sex; 
Sex_all2    = [Sex_HCs; Sex_disPTs]; 
Sex_all2_num = double(strcmp(Sex_all2,'M'));
Sex_all4     = [Sex_HCs; Sex_disPTs; Sex_valPTs; Sex_sesPTs]; 
Sex_all4_num = double(strcmp(Sex_all4,'M'));

Side_allPTs = [Demo_disPTs.Side; Demo_valPTs.Side; Demo_senPTs.Side];

Dur_PTs = Demo_disPTs.Duration;

SZFreq_PTs = Demo_disPTs.SZFreq;   
SZFreq_PTs_num = zeros(size(SZFreq_PTs));
SZFreq_PTs_num(strcmp(SZFreq_PTs, 'Daily'))   = 1;   
SZFreq_PTs_num(strcmp(SZFreq_PTs, 'Weekly'))  = 2;
SZFreq_PTs_num(strcmp(SZFreq_PTs, 'Monthly')) = 3; 
SZFreq_PTs_num(strcmp(SZFreq_PTs, 'Yearly'))  = 4;

Prog_PTs = Demo_disPTs.x1YProgI; 
Prog_PTs_bi = repmat({'NSF'}, size(Prog_PTs)); 
Prog_PTs_bi(strcmp(Prog_PTs, 'ILAE 1')) = {'SF'};
% Prog_num = double(strcmp(Prog_PTs, 'ILAE 1'));