function [Qnew,Snew,SpillNew] = CheckStorage(Q,Evdp,I,VarMax,VarMin)
  nVar = numel(Q);
  chckArray = zeros(1,nVar);
  [S,Spill] = StorageEval(Q,Evdp,I);
  
  while chckArray == S
     
      for t = 1: nVar
            %Q(t) = unifrnd(VarMin(t), VarMax(t));
                Q(t)=VarMin(t)+(VarMax(t)-VarMin(t)).*rand(size(VarMin(t)));
      end
        [S,Spill] = StorageEval(Q,Evdp,I);

  end
  
  Qnew = Q;
  Snew = reshape(S,1,nVar);
  SpillNew = reshape(Spill,1,nVar);

end