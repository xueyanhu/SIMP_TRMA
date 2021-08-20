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
%%% This is a code for  SIMP with its optimization formulation:
%%%				find: rou=[1e-3,1]
%%%				min:   F*U
%%%				s.t.
%%%					V-Vf <= 0;
%%% Fliter Type: same as 88Lines_SIMP
%%%
%%% CASE: Canti Beam
%%% METHOD: TRMA_Mixed
%%%==========================================================================%%%
clc
clear all
%%% <<Function Lab>>
% 1: Laptop %2: AAU workstation %3: ZJU workstation
addpath(genpath(pwd));

%%% <<Initialization>>
%% - SIMP par
global pSIMP; pSIMP = 3;
%% - Mesh Par (The default element side length is 1)
for nelx = [100]
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
    vf = 0.5, type1 = 'ineq';
    nCons1 = 1;
    consInfo1.val = 0.0; consInfo1.sen = zeros(1,nCons1);
    consInfo1.lam = 1000*ones(nCons1,1); consInfo1.type = type1;
    consInfo1.cof = 1;
    ALConsInfo = {};
    consInfo{1} = consInfo1;
    % - Initialize penality value and control parameter
    sig = 1100
    lamMax = 90000; sigMax = 10000;
    % - Initialize radius of TR and limits of design variables
    xmax = 1.0*ones(nDV,1); xmin = 0*ones(nDV,1);
    dBound = 0.02*ones(nDV,1);
    dBoundMax = min(abs(xmax-rou),abs(rou-xmin));

    %% - TRMA parameter Initialization
    global TRMAPar

    TRMAPar.criGama1 = 0.1; TRMAPar.criGama2 = 0.9;
    TRMAPar.fliTol = 1e-4;
    %% - Convergence Parameters
    Tol_d = 1e-4;
    Tol_h = 1e-4;
    iter = 0; iter_h = 0; iter_d = 0;
    %% -Initialize Records
    objRecord = []; consRecord = []; DVsR = [rou];
    rouViola = [];

    %%% <<First Loop>>
    [Eele,rouFli,dEele] = Map_Line88(rou,H,Hs,materialInfo);
    FEAResult = PDESolver_FEA_MATLAB(meshInfo,Eele,materialInfo,...
        BCInfo,preMatrixInfo,OutputRequest);
    objInfo = OBJ_SimpComplianceLine88(FEAResult,dEele,H,Hs);
    consInfo = CONS_SimpVolumeLine88(rouFli,vf,nelx,nely,consInfo,H,Hs);
    %- Plot & Print State
    colormap(gray); imagesc(1-reshape(rou,nely,nelx)); caxis([0,1]); 
    axis equal, axis off;pause(1e-6);
    [d,modelInfo] = TRMA_SolSub(nDV,objInfo,consInfo,ALConsInfo,...
        sig,dBound,xmax,xmin);
    modelInfo1=modelInfo;
    dIndi = sqrt(d'*d);
    maxCons = TRMA_calMaxCons(consInfo,ALConsInfo);
    lamCri = TRMA_calMaxLam(consInfo,ALConsInfo);
    statePrint(iter,iter_h,iter_d,objInfo,maxCons,dIndi,lamCri,sig);
    h = TRMA_calH(consInfo,ALConsInfo);
    h=1;
    hRecord = [h];
    Filter = [h,objInfo.val; 10*h,inf];
    TRMAPar.outC1 = 1;TRMAPar.outC2 = 0.05;
    TRMAPar.inC1 = 1;TRMAPar.inC2 = 0.05;;
    while  h > Tol_h   && iter < 2000;
        iter_h = iter_h+1;
        dIndi = sqrt(d'*d);
        iter_d = 0; lamCri = TRMA_calMaxLam(consInfo,ALConsInfo);
        while dIndi > Tol_d   && sig<sigMax && lamCri<lamMax && iter < 2000
            iter_d = iter_d+1;iter = iter+1;
            % - Giving a trial step
            rouTrial = rou+d;
            % - Update the state at trial step
            [EeleTrial,rouFliTrial,dEeleTrial] = ...
                Map_Line88(rouTrial,H,Hs,materialInfo);
            FEAResultTrial = PDESolver_FEA_MATLAB(meshInfo,EeleTrial,materialInfo...
                ,BCInfo,preMatrixInfo,OutputRequest);
            objInfoTrial = OBJ_SimpComplianceLine88(FEAResultTrial,dEeleTrial,H,Hs);
            consInfoTrial = CONS_SimpVolumeLine88(rouFliTrial,vf,nelx,nely,...
                consInfo,H,Hs);
            ALConsInfoTrial = {};
            % - Update the lam at trial step
            [consInfoTrial,ALConsInfoTrial] = TRMA_calLamTrial(sig,consInfoTrial,...
                ALConsInfoTrial);
            % - Calculate criterion
            AL0 = TRMA_GenAL(objInfo,consInfo,ALConsInfo,sig);
            ALTrial = TRMA_GenAL(objInfoTrial,consInfoTrial,ALConsInfoTrial,sig);
            AlMol0  = AL0;
            AlMolD = modelInfo1.val;
            % - Judge
            [AccepCri,cri] = TRMA_AcceByCri(ALTrial,AL0,AlMol0,AlMolD);
            % - Relaxition Judge by Fliter
            [AccepFli,Filter] = TRMA_AcceByFli(objInfoTrial,consInfoTrial,...
                ALConsInfoTrial,Filter);
            Accep = AccepCri | AccepFli;

            % - Update state by criterion
            dBoundMax = min(abs(xmax-rouTrial),abs(rouTrial-xmin));
            [rou,consInfo,ALConsInfo,objInfo,sig,dBound] = ...
            TRMA_UpDate(Accep,cri,objInfo,objInfoTrial,consInfo,consInfoTrial,...
                ALConsInfo,ALConsInfoTrial,sig,dBound,d,dBoundMax,rou,xmax,xmin);
            if Accep
                rouFli = rouFliTrial; Eele = EeleTrial; dEele = dEeleTrial;
                FEAResult = FEAResultTrial;
            else
                rouFli = rouFli; Eele = Eele; dEele = dEele;
                FEAResult = FEAResult;
            end
            [d,modelInfo] = TRMA_SolSub(nDV,objInfo,consInfo,ALConsInfo,...
                sig,dBound,xmax,xmin);
            modelInfo1=modelInfo;
            dIndi = sqrt(d'*d);
            %- Plot & Print State
            rou = max(round(rou*1e4)/1e4,0);
            rouFli = max(round(rouFli*1e4)/1e4,0);
            maxCons = TRMA_calMaxCons(consInfo,ALConsInfo);
            lamCri = TRMA_calMaxLam(consInfo,ALConsInfo);
            colormap(gray); imagesc(1-reshape(rou,nely,nelx)); caxis([0,1]);
            axis equal, axis off; pause(1e-6);
            statePrint(iter,iter_h,iter_d,objInfo,maxCons,dIndi,lamCri,sig)
            %- Record state
            objRecord = [objRecord;objInfo.val];
            consRecord = [consRecord;maxCons];
            denSta = (sum(rou>=1-1e-4)+sum(rou<=1e-3+1e-4))/length(rou);
            disp(['Den:' sprintf('%6.3f\t',1-denSta)])
            rouViola = [rouViola,denSta];
            h = TRMA_calH(consInfo,ALConsInfo);
            disp(['h:' sprintf('%6.3f\t',h)])
        end
        h = TRMA_calH(consInfo,ALConsInfo);
        if h > TRMAPar.outC1
            sig = 1.5*sig;
        elseif  h < TRMAPar.outC2
            %TRMAPar.outC1 = max(TRMAPar.outC1/1.2, 2*TRMAPar.outC2);
            TRMAPar.outC2= max(min(TRMAPar.inC2/1.2,0.08*h),1e-5);
            TRMAPar.inC2= max(min(TRMAPar.inC2/1.2,0.08*h),1e-5);
            sig = sig/2;
        else
            sig = 1.05*sig;
        end
        consInfo1 = consInfo{1};
        lam = 1*consInfo1.lam;
        [consInfo,ALConsInfo] = TRMA_calLamTrial(sig,consInfo,...
                ALConsInfo);
        consInfo1 = consInfo{1};
        consInfo1.lam = lam+0.6*(consInfo1.lam-lam);
        consInfo{1} = consInfo1;
        lamMax = max(10*consInfo1.lam,10);
        [Eele,rouFli,dEele] = Map_Line88(rou,H,Hs,materialInfo);
        FEAResult = PDESolver_FEA_MATLAB(meshInfo,Eele,materialInfo,...
            BCInfo,preMatrixInfo,OutputRequest);
        objInfo = OBJ_SimpComplianceLine88(FEAResult,dEele,H,Hs);
        consInfo = CONS_SimpVolumeLine88(rouFli,vf,nelx,nely,consInfo,H,Hs);
        dBound = max(dBoundMax/2,1e-6);
        [d,modelInfo] = TRMA_SolSub(nDV,objInfo,consInfo,ALConsInfo,...
        sig,dBound,xmax,xmin);
        modelInfo1 = modelInfo;
        dIndi = sqrt(d'*d);
        sigMax = 10*sig;
        iter = iter+1;
        colormap(gray); imagesc(1-reshape(rou,nely,nelx)); caxis([0,1]);
            axis equal, axis off; pause(1e-6);
        statePrint(iter,iter_h,iter_d,objInfo,maxCons,dIndi,lamCri,sig)
    end
    fileName = ['SIMPTRMA',num2str(nelx)]
    save(fileName)
end
