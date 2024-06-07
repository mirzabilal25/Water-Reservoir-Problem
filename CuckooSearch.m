
function out = CuckooSearch(problem, CSparams)

%% Problem Definition
CostFunction = problem.CostFunction;
 
 nVar = problem.nVar;
 VarSize = problem.VarSize;
 VarMin = problem.VarMin;
 VarMax = problem.VarMax;
 
Demand = problem.Demand;
Evdp = problem.Evdp;
Inflow = problem.Inflow;
FEMax = problem.FEMax;
FE = 0;
FECost = [];
BestCostVal = inf;
 
%% Parameter Declaration

% Number of nests (or different solutions)
n=CSparams.n;
% Discovery rate of alien eggs/solutions
pa=CSparams.pa;
% Tolerance
Tol=CSparams.Tol;

%nd=15;
nd =CSparams.nd;
% Lower bounds 
Lb=CSparams.Lb;     %-5*ones(1,nd);
%Lb=-10;
% Upper bounds
Ub=CSparams.Ub;    %5*ones(1,nd);
%Ub = 10;

MaxIt =CSparams.MaxIt;
FE = 0;
FECost = [];
BestCostVal = inf;
%% Initialization

% Random initial solutions
for i=1:n,
    %nest(i,:)=Lb+(Ub-Lb).*rand(size(Lb));
    %nest(i,:) = unifrnd(VarMin, VarMax,VarSize);
    for t = 1: nVar
       nest(i,t) = VarMin(t)+(VarMax(t)-VarMin(t)).*rand(size(VarMin(t)));
    end
end

% Get the current best
%fitness=10^10*ones(n,1);
fitness=10^10*ones(n,1);
newnest=nest;
for j=1:size(nest,1),
    [newnest(j,:),Storage(j,:),Spill(j,:)] = CheckStorage(newnest(j,:),Evdp,Inflow,VarMax,VarMin);
    [fnew,FE]=CostFunction(newnest(j,:),Demand,FE);
    if (FE == 1)
        FECost(FE) = (fnew);
    else
        FECost(FE) = (fnew < BestCostVal)*fnew + (fnew >= BestCostVal)*BestCostVal;
    end
    BestCostVal = FECost(FE);
    
    if fnew<=fitness(j),
       fitness(j)=fnew;
       nest(j,:)=newnest(j,:);
    end
end
% Find the current best
[fmin,K]=min(fitness) ;
best.Position=nest(K,:);
best.Storage = Storage(K,:);
best.Spill = Spill(K,:);
best.Cost = fmin;
bestnest = best;


BestCosts = zeros(MaxIt,1);
%% Main Loop Cuckoo Search 
it = 1;
BCost = inf;
%for iter=1:MaxIt
%while it<=MaxIt && BCost > 0
while FE <= FEMax

    % Generate new solutions (but keep the current best)
     %new_nest=get_cuckoos(nest,bestnest,Lb,Ub);   
     n=size(nest,1);
     beta=3/2;
     sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
     for j=1:n,
        s=nest(j,:);
        u=randn(size(s))*sigma;
        v=randn(size(s));
        step=u./abs(v).^(1/beta);
        stepsize=0.01*step.*(s-bestnest.Position);
        s=s+stepsize.*randn(size(s));
        new_nest(j,:)=simplebounds(s,Lb,Ub);
     end
     [fnew,best,nest,fitness,FE,BestCostVal,FECost]=get_best_nest(nest,new_nest,fitness,CostFunction,Demand,Evdp,Inflow,VarMax,VarMin,FE,BestCostVal,FECost);
     % Update the counter
      %N_iter=N_iter+n; 
      
    % Discovery and randomization
      %new_nest=empty_nests(nest,Lb,Ub,pa) ;
      n=size(nest,1);
      K=rand(size(nest))>pa;
        stepsize=rand*(nest(randperm(n),:)-nest(randperm(n),:));
        new_nest=nest+stepsize.*K;
        for j=1:size(new_nest,1)
            s=new_nest(j,:);
            new_nest(j,:)=simplebounds(s,Lb,Ub);  
        end
    % Evaluate this set of solutions
      [fnew,best,nest,fitness,FE,BestCostVal,FECost]=get_best_nest(nest,new_nest,fitness,CostFunction,Demand,Evdp,Inflow,VarMax,VarMin,FE,BestCostVal,FECost);
    % Update the counter again
      %N_iter=N_iter+n;
      
    % Find the best objective so far  
    if fnew<fmin,
        fmin=fnew;
        bestnest=best;
    end
    BestCosts(it) = fmin;
    BCost = fmin;
    %disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]);
    it = it+1;
end
    disp(['Total Iteration: ' num2str(it-1)]);
    out.BestSol = bestnest;
    out.BestCost = fmin;
    out.BestCostsCS = BestCosts;
    out.FE = FE;
    out.FECost = FECost;
    out.iter = it-1;
end

%{
%% Display all the nests
%disp(strcat('Total number of iterations=',num2str(N_iter)));
%fmin
%bestnest

figure,
semilogy(BestCosts, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost CS');
grid on;

disp(['Best Sol Position:' num2str(fmin)]);
disp(['Solution Cost:' num2str(bestnest)]);
%}

