##  cyclic shifts

Read("coxeter.g");

cyclic_shifts:= function(W, w)
    local   byCyclicShift;
    byCyclicShift:= function(x, s)
        local   lx,  y,  ly;
        lx:= coxeterLength(W, x);
        y:= x^s;#
        ly:= coxeterLength(W, y);
        if lx = ly then  return y;  fi;
        return x;
    end;
    return orbit(GeneratorsOfGroup(W), w, byCyclicShift);
end;

##  the cyclic shift class of an involution x is just {x}

cyclic_shifts_with_edges:= function(W, w)
    local   byCyclicShift;
    byCyclicShift:= function(x, s)
        local   lx,  y,  ly;
        lx:= coxeterLength(W, x);
        y:= x^s;#
        ly:= coxeterLength(W, y);
        if lx = ly then  return y;  fi;
        return x;
    end;
    return orbit_with_edges(GeneratorsOfGroup(W), w, byCyclicShift);
end;

