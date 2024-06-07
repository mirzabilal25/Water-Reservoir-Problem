%% Find the current best nest
function [fmin,bestN,nest,fitness,FEnew,BestCostVal,FECost]=get_best_nest(nest,newnest,fitness,CostFunction,Demand,Evdp,Inflow,VarMax,VarMin,FE,BestCostVal,FECost)
    FEnew = FE;
    % Evaluating all new solutions
for j=1:size(nest,1),
    [newnest(j,:),Storage(j,:),Spill(j,:)] = CheckStorage(newnest(j,:),Evdp,Inflow,VarMax,VarMin);
    [fnew,FEnew]=CostFunction(newnest(j,:),Demand,FEnew);
    FECost(FEnew) =  (fnew < BestCostVal)*fnew + (fnew >= BestCostVal)*BestCostVal;
    BestCostVal = FECost(FEnew);
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
bestN = best;

