onInvolutions := function(x, s)
    if x * s = s * x then return OnRight(x, s); fi;
    return OnPoints(x, s);
end;

onInvolutionClasses:= function(x, a)
    local y;
    y:= Representative(x);
    if y^a <> y then return x; fi;
    return OnRight(y, a)^ActingDomain(x);
end;

involutionClasses := function(group)
    local gens;
    gens:= GeneratorsOfGroup(group);
    return orbit(orbits(gens, gens, OnPoints), Identity(group)^group, onInvolutionClasses);
end;
