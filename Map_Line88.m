function [Eele,rouFli,dEele] = Map_Line88(DV,H,HS,materialInfo);
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
%%% This is a Map function to fliter the element density and map it into the
%%% element Modulus by the SIMP model shown in the 88-line code
%%%
%%%==========================================================================%%%
global pSIMP

Es = materialInfo.Es; Emin = materialInfo.Emin;
rouFli = (H*DV)./HS;
Eele = Emin+rouFli.^pSIMP*(Es-Emin);
dEele = pSIMP*rouFli.^(pSIMP-1)*(Es-Emin);
