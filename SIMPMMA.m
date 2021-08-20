%%%===========================Copyright======================================%%%
%%%   Version August. 2021
%%%
%%%   Xueyan Hu <huxueyan@zju.edu.cn>
%%%   PHD student in
%%%   	Institute of Applied Mechanics,Zhejiang University
%%%	  Guest Student in
%%%   	Department of Mechanical and Manufacturing Engineering ,
%%%											Aalborg University
%%%
%%%===========================Description====================================%%%
%%% This is a code for  PISIMP with its optimization formulation:
%%%				find: rou=[1e-3,1]
%%%				min:   F*U
%%%				s.t.
%%%					V-Vf <= 0;
%%%                 rou*(1-rou)=0
%%% Fliter Type: same as 88Lines_SIMP
%%%
%%% CASE: Canti Beam
%%% METHOD: TRMA_Mixed
%%%==========================================================================%%%
clc
clear all
%%% <<Function Lab>>
% 1: Laptop %2: AAU workstation %3: ZJU workstation
computerNumber = 1;
if computerNumber == 1;
	addpath(genpath(...
		'E:\OneDrive - zju.edu.cn\code\TopOptFunLab_Ver3.0'));
elseif computerNumber == 2;
	addpath(genpath(...
		'D:\xuhu\code0727\TopOptFunLab_Ver3.0'));
elseif computerNumber == 3;
	addpath(genpath(...
		'E:\xuhu\code\code\TopOptFunLab_Ver3.0'));
end
%%% <<Initialization>>
%% - SIMP par
global pSIMP; pSIMP = 3;
%% - Mesh Par (The default element side length is 1)
for nelx = [100]%,600,800]
    nely = nelx/2;
    l = nelx; h = nely;
    elementType = 'CPS4'
    %% - Material par
    materialInfo.Es = 1.0; materialInfo.Emin = 1e-9;
    materialInfo.nu = 0.3;

    %%- MESH
    [meshInfo,contourShow] = Init_PDEslover_Mesh_UniforRec_CPS4_1Size(nelx,nely,...
        elementType);

    %% - BC condition
    coord = meshInfo.coord;
    % - fix Dof
    bound = find(coord(:,1)==0); fixdofs = sort([2*bound-1;2*bound]);
    % -LOAD
    F = sparse(2*(nely+1)*(nelx+1),1); Fc=-1;
    loadPoint = find(coord(:,1)==l&coord(:,2)==h/2); F(2*loadPoint,1) = Fc;
    BCInfo.fixdof = fixdofs; BCInfo.F = F;

    %% - Output request
    OutputRequest.History{1} = 'Comp';
    OutputRequest.FieldAtEle{1} = 'EleCompSolid';

    %% - pre-process
    [preMatrixInfo] = PDE_preProcess_matrixInfo(meshInfo,elementType,materialInfo,...
        OutputRequest);
    %% - Initialize  Design variable
    nDV = meshInfo.nEl;
    rou = 0.501*ones(nDV,1);
    %%- Fliter Initialization
    rmin=1.5;
    [H,Hs] = Filter_line88(nelx,nely,rmin);
    %%- Initialize constriant format
    % - volume constraint
    vf = 0.5;
    nCons1 = 1;
    consInfo1.val = 0.0; consInfo1.sen = zeros(nDV,nCons1);
    consInfo{1} = consInfo1;
    ALConsInfo = {};
    % - Initialize limits of design variables
    xmax = 1.0*ones(nDV,1); xmin = 0*ones(nDV,1);

    % - Initialize MMA
    xval=rou; xold1 = rou; xold2 = rou;
    low   = xmin;
	upp   = xmax;
    c=1000*ones(nCons1,1); d=ones(nCons1,1);
	a0=1; a=zeros(nCons1,1);
	moLi = 0.8;

	% - Initialize Records
	objRecord = []; consRecord = []; rouViola = [];

	% - Initialize Tol
	change = 1; iter = 0; tol = 1e-4; maxIter = 2000;objOld = 1000;

	%%% <<Main Loop>>
	while change > tol && iter <= maxIter
		[Eele,rouFli,dEele] = Map_Line88(rou,H,Hs,materialInfo);
		FEAResult = PDESolver_FEA_MATLAB(meshInfo,Eele,materialInfo,...
        BCInfo,preMatrixInfo,OutputRequest);
        objInfo = OBJ_SimpComplianceLine88(FEAResult,dEele,H,Hs);
    	consInfo = CONS_SimpVolumeLine88(rouFli,vf,nelx,nely,consInfo,H,Hs);
    	%- Plot & Print State
	    colormap(gray); imagesc(1-reshape(rou,nely,nelx)); caxis([0,1]); 
	    axis equal, axis off;pause(1e-6);
	    maxCons = TRMA_calMaxCons(consInfo,ALConsInfo);
	    %- Record state
        objRecord = [objRecord;objInfo.val];
        consRecord = [consRecord;maxCons];
        denSta = (sum(rou>=1-1e-4)+sum(rou<=0+1e-4))/length(rou);
        rouViola = [rouViola,denSta];

        %- MMA process
        f0val = objInfo.val; df0dx = objInfo.sen; %df0dx = df0dx/max(abs(df0dx));
        consInfo1 = consInfo{1};
        fval(1,1) = consInfo1.val; dfdx = (consInfo1.sen)'; %dfdx = dfdx/max(abs(dfdx));

        %- Optimizer MMA_Version Fixed Move Limit
        [xmma,ymma,zmma,lam,xsi,eta,mu,zet,ss,low,upp] = ...
        OptSolver_mmasubFixMoLi(nCons1,nDV,iter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d,moLi);

        % - Update state
        iter = iter+1;
        xold2 = xold1; xold1 = xval; xval = xmma; rou = (xval*1e4)/1e4;
        change = abs(f0val-objOld)/abs(objOld);
        objOld = f0val;
        %--State print
   		 disp([' It.: ' sprintf('%4i\t',iter) ' Obj.: ' sprintf('%6.3f\t',f0val) ' Vol.: ' ...
        sprintf('%6.4f\t',fval) 'ch.:' sprintf('%6.4f\t',change)]);
   	end
 end




