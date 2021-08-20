function  objInfo = OBJ_SimpComplianceLine88(FEAResult,dEele,H,HS)
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
%%% This function is to generate the information of objective function including
%%% value at current step 'objInfo.val' and the sensitivity at current step
%%% 'objInfo.sen'
%%%
%%%==========================================================================%%%
global pSIMP;
eleCompSolid = FEAResult.FieldAtEle.EleCompSolid;
c = FEAResult.History.Comp;
dc = -dEele.*eleCompSolid;
dc = H*(dc./HS);

objInfo.val = c; objInfo.sen = dc;