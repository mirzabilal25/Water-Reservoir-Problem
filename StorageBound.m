function [Snew,Spill] = StorageBound(Sval)
  Snew = Sval;
  Spill = 0;
    if Sval < 0
        Snew = inf;
    elseif Sval > 608
        Spill = Sval - 608;
        Snew = 608;
    end

end