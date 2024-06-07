clc
clear;
close all;

%% Problem Definition
 MaxRun = 2;   % Number of runs
 MaxIter =500; % Number of iterations
 popSize = 100; % population size
 problem.CostFunction = @(Q,Demand,FE) waterReleaseOpFunc(Q,Demand,FE);
 
[problem.Demand,problem.Evdp,problem.Inflow] = initProblemVariabls(); 
 problem.nVar = 12;  % number of decision variables, just for 1 year
 problem.VarSize=[1 problem.nVar];
 problem.VarMin = zeros(1,problem.nVar);
 problem.VarMax = problem.Demand; 
 problem.Storage = zeros(1,problem.nVar);
 problem.FEMax = 5000;          % maximum number of function evaluation
 FE = [];
%% PSO Optimization

% Parameters of PSO
PSOparams.MaxIt = MaxIter;
PSOparams.nPop =  popSize;
PSOparams.w=1; % weight
PSOparams.wdamp = 0.99; % adaptive weight
PSOparams.c1 =1.494;    
PSOparams.c2 =1.494; 

% Calling PSO
for r = 1: MaxRun
    tic
    outPSO = PSO(problem,PSOparams);
    timePSO(r) = toc;
    BestSolnPSO(r) = outPSO.BestSolution;
    PSOCost(r) = BestSolnPSO(r).Cost;
    FE(1) = outPSO.FE;
    BestCostsPSO(r,:) = outPSO.BestCostsPSO;
    FECostPSO(r,:) = outPSO.FECost;
    iter(1) = outPSO.iter;
end
[Cost,index]= min(PSOCost);
OptCost(1) = Cost;
PSOBest = BestSolnPSO(index);
PSOBestCosts = BestCostsPSO(index,:);
PSOBestCostsFE = FECostPSO(index,:);
PSOtimeEstd =  timePSO(index);

disp(['FE PSO:' num2str(outPSO.FE)]);
disp(['Best Sol Position PSO:' num2str(PSOBest.Position)]);
disp(['Best Sol Storage PSO:' num2str(PSOBest.Storage)]);
disp(['Best Sol Spill PSO:' num2str(PSOBest.Spill)]);
disp(['Solution Cost PSO:' num2str(PSOBest.Cost)]);
disp(['Estimated Time: ' num2str(PSOtimeEstd)]);
disp(char(10));
%
%% DE Optimization

% Parameters of DE
DEparams.MaxIt=MaxIter;      % Maximum Number of Iterations
DEparams.nPop=popSize;            % Population Size
DEparams.beta_min=0.2;       % Lower Bound of Scaling Factor
DEparams.beta_max=0.8;       % Upper Bound of Scaling Factor
DEparams.pCR=0.2;

% Call DE
for r = 1: MaxRun
    tic
    outDE = DEImp(problem, DEparams);
    timeDE(r) = toc;
    BestSolDE(r) = outDE.BestSol;
    DECost(r) = BestSolDE(r).Cost;
    FE(2) = outDE.FE;
    BestCostsDE(r,:) = outDE.BestCostsDE;
    FECostDE(r,:) = outDE.FECost;
    iter(2) = outDE.iter;
end
[Cost,index]= min(DECost);
OptCost(2) = Cost;
DEBest = BestSolDE(index);
DEBestCosts = BestCostsDE(index,:);
DEBestCostsFE = FECostDE(index,:);
DEtimeEstd =  timeDE(index);

disp(['FE DE:' num2str(outDE.FE)]);
disp(['Best Sol Position DE:' num2str(DEBest.Position)]);
disp(['Best Sol Storage DE:' num2str(DEBest.Storage)]);
disp(['Best Sol Spill DE:' num2str(DEBest.Spill)]);
disp(['Solution Cost DE:' num2str(DEBest.Cost)]);
disp(['Estimated Time: ' num2str(DEtimeEstd)]);
disp(char(10));

%% ABC Optimization

% Parameters of ABC 
ABCparams.MaxIt=MaxIter;                                            % Maximum Number of Iterations
ABCparams.nPop=popSize;                                                  % Population Size (Colony Size)
ABCparams.nOnlooker=ABCparams.nPop;                                 % Number of Onlooker Bees
ABCparams.L=round(0.6*problem.nVar*ABCparams.nPop);                 % Abandonment Limit Parameter (Trial Limit)
ABCparams.a=1;                                                      % Acceleration Coefficient Upper Bound

% Call ABC
for r = 1: MaxRun
    tic
    outABC = ABC(problem, ABCparams);
    timeABC(r) = toc;
    BestSolABC(r) = outABC.BestSol;
    ABCCost(r) = BestSolABC(r).Cost;
    FE(3) = outABC.FE;
    BestCostsABC(r,:) = outABC.BestCostsABC;
    FECostABC(r,:) = outABC.FECost;
    iter(3) = outABC.iter;
end
[Cost,index]= min(ABCCost);
OptCost(3) = Cost;
ABCBest = BestSolABC(index);
ABCBestCosts = BestCostsABC(index,:);
ABCBestCostsFE = FECostABC(index,:);
ABCtimeEstd =  timeABC(index);


disp(['FE ABC:' num2str(outABC.FE)]);
disp(['Best Sol Position ABC:' num2str(ABCBest.Position)]);
disp(['Best Sol Storage ABC:' num2str(ABCBest.Storage)]);
disp(['Best Sol Spill ABC:' num2str(ABCBest.Spill)]);
disp(['Solution Cost ABC:' num2str(ABCBest.Cost)]);
disp(['Estimated Time: ' num2str(ABCtimeEstd)]);
disp(char(10));
%% GA Optimization

% GA Parameters
GAparams.MaxIt=MaxIter;                                 % Maximum Number of Iterations
GAparams.nPop=popSize;                                       % Population Size
GAparams.pc=0.7;                                        % Crossover Percentage
GAparams.nc=2*round(GAparams.pc*GAparams.nPop/2);       % Number of Offsprings (also Parnets)
GAparams.gamma=0.4;                                     % Extra Range Factor for Crossover
GAparams.pm=0.3;                                        % Mutation Percentage
GAparams.nm=round(GAparams.pm*GAparams.nPop);           % Number of Mutants
GAparams.mu=0.1;                                        % Mutation Rate
% RouletteWheelSelection
GAparams.beta=8;                                        % Selection Pressure

% Call GA
for r = 1: MaxRun
    tic
    outGA = GA(problem, GAparams);
    timeGA(r) = toc;
    BestSolGA(r) = outGA.BestSol;
    GACost(r) = BestSolGA(r).Cost;
    FE(4) = outGA.FE;
    BestCostsGA(r,:) = outGA.BestCostsGA;
    FECostGA(r,:) = outGA.FECost;
    iter(4) = outGA.iter;
end
[Cost,index]= min(GACost);
OptCost(4) = Cost;
GABest = BestSolGA(index);
GABestCosts = BestCostsGA(index,:);
GABestCostsFE = FECostGA(index,:);
GAtimeEstd =  timeGA(index);

disp(['FE GA:' num2str(outGA.FE)]);
disp(['Best Sol Position GA:' num2str(GABest.Position)]);
disp(['Best Sol Storage GA:' num2str(GABest.Storage)]);
disp(['Best Sol Spill GA:' num2str(GABest.Spill)]);
disp(['Solution Cost GA:' num2str(GABest.Cost)]);
disp(char(10));

%% Cuckoo Search Optimization

% Parameter Declaration
CSparams.n=popSize;                        % Number of nests (or different solutions)
CSparams.pa=0.25;                     % Discovery rate of alien eggs/solutions
CSparams.Tol=1.0e-5;                  % Tolerance
CSparams.nd = 12;
CSparams.Lb=zeros(1,problem.nVar);    %-5*ones(1,nd); % Lower bounds
CSparams.Ub=problem.Demand;    % Upper bounds
CSparams.MaxIt = MaxIter;

% Call CS
for r = 1: MaxRun
    tic
    outCS = CuckooSearch(problem, CSparams);
    timeCS(r) = toc;
    BestCostCS(r) = outCS.BestCost;
    BestSolCS(r) = outCS.BestSol;
    FE(5) = outCS.FE;
    BestCostsCS(r,:) = outCS.BestCostsCS;
    FECostCS(r,:) = outCS.FECost;
    iter(5) = outCS.iter;
end
[CSCost,index]= min(BestCostCS);
OptCost(5) = CSCost;
CSBest = BestSolCS(index);
CSBestCosts = BestCostsCS(index,:);
CSBestCostsFE = FECostCS(index,:);
CStimeEstd =  timeCS(index);

disp(['FE CS:' num2str(outCS.FE)]);
disp(['Best Sol Position CS:' num2str(CSBest.Position)]);
disp(['Best Sol Storage CS:' num2str(CSBest.Storage)]);
disp(['Best Sol Spill CS:' num2str(CSBest.Spill)]);
disp(['Solution Cost CS:' num2str(CSCost)]);
disp(['Estimated Time: ' num2str(CStimeEstd)]);
disp(char(10));

%% SA Optimization

% Parameters of SA

SAparams.MaxIt=MaxIter;        % Maximum Number of Iterations
SAparams.MaxSubIt=65;       % Maximum Number of Sub-iterations
SAparams.T0=1; %0.025;          % Initial Temp.
SAparams.alpha=0.8;
SAparams.Tmin = 1e-70;

% Call SA
for r = 1: MaxRun
    tic
    outSA = SA(problem, SAparams);
    timeSA(r) = toc;
    BestSolSA(r) = outSA.BestSol;
    SACost(r) = BestSolSA(r).Cost;
    FE(6) = outSA.FE;
    BestCostsSA(r,:) = outSA.BestCostsSA;
    FECostSA(r,:) = outSA.FECost;
    iter(6) = outSA.iter;
end
[Cost,index]= min(SACost);
OptCost(6) = Cost;
SABest = BestSolSA(index);
SABestCosts = BestCostsSA(index,:);
SABestCostsFE = FECostSA(index,:);
SAtimeEstd =  timeSA(index);

disp(['FE SA:' num2str(outSA.FE)]);
disp(['Best Sol Position SA:' num2str(SABest.Position)]);
disp(['Best Sol Storage SA:' num2str(SABest.Storage)]);
disp(['Best Sol Spill SA:' num2str(SABest.Spill)]);
disp(['Solution Cost SA:' num2str(SABest.Cost)]);
disp(['Estimated Time: ' num2str(SAtimeEstd)]);
disp(char(10));

%% RESULT
MaxIter = max(iter);
MaxFuncEval = max(FE);

figure;
plot(PSOBestCostsFE);
hold on
plot(DEBestCostsFE);
hold on
plot(ABCBestCostsFE);
hold on
plot(GABestCostsFE);
hold on
plot(CSBestCostsFE);
hold on
plot(SABestCostsFE);
axis([0 MaxFuncEval 0 inf]);
legend('PSO','DE','ABC','GA','CS','SA');
hold off
xlabel('Function Evaluation');
ylabel('Best Cost');
grid on;

%{
figure,
bar(FE);
set(gca,'xticklabel',{'PSO','DE','ABC','GA','CS','SA'});
%}

%
% For Iteration
figure;
plot(PSOBestCosts);
hold on
plot(DEBestCosts);
hold on
plot(ABCBestCosts);
hold on
plot(GABestCosts);
hold on
plot(CSBestCosts);
hold on
plot(SABestCosts);
axis([0 MaxIter 0 inf]);
legend('PSO','DE','ABC','GA','CS','SA');
hold off
xlabel('Iteration');
ylabel('Best Cost');
grid on;
%}