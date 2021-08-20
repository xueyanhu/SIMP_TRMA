function consInfoValSenUp = CONS_SimpVolumeLine88(rouFli,vf,nelx,nely,consInfo,H,HS);
%%%===========================Copyright======================================%%%
%%%   Version July. 2021
%%%
%%%   Xueyan Hu <huxueyan@zju.edu.cn>
%%%   PHD student in
%%%   	Institute of Applied Mechanics,Zhejiang University
%%%	  Guest Student in
%%%   	Department of Mechanical and Manufacturing Engineering ,
%%%											Aalborg University
%%%
%%%===========================Description====================================%%%
%%% This function is to generate the information of constriant  including
%%% value at current step 'consInfo{i}.val' and the sensitivity at current step
%%% 'consInfo{i}.sen'
%%%
%%%==========================================================================%%%
%% - volume constraint
Vval = mean(rouFli)-vf;%*nelx*nely;
Vsen = ones(size(rouFli))/length(rouFli);
Vsen = H*(Vsen./HS);
consInfo1 = consInfo{1};
consInfo1.val = Vval; consInfo1.sen = Vsen; consInfoValSenUp{1} = consInfo1;