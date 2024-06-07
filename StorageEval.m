function [S,Spill] = StorageEval(Q,Evdp,I)


nVar = numel(Q);
evap = zeros(1, nVar);
S = zeros(1, nVar);
Spill = zeros(1, nVar);

Q = reshape(Q,1,nVar);
I = reshape(I,1,nVar);
Evdp = reshape(Evdp,1,nVar);
%% Main processing

for t = 2: nVar
    S(t) = S(t-1)+I(t-1)-Q(t-1)-evap(t-1);
    [S(t),Spill(t)] = StorageBound(S(t));
    if S(t) == inf
        S = zeros(1,nVar);
        break;
    end
    av_st  = S(t) + S(t-1)/2;
    surf_area = 16.025 + 0.0854*av_st - 4*(10^-5)*av_st^2 + (10^-8)*(av_st)^3;
    evap(t) = Evdp(t-1).*surf_area;
end
  


end