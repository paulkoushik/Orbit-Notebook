onInvolutions := function(x, s)
    if x * s = s * x then return OnRight(x, s); fi;
    return OnPoints(x, s);
end;

involutions := function(W)
    return orbit(GeneratorsOfGroup(W), Identity(W), onInvolutions);
end;

onInvolutionClasses := function(x, a)
    local y;
    y := Representative(x);
    if y^a <> y then return x; fi;
    return OnRight(y, a)^ActingDomain(x);
end;

involutionClasses := function(W)
    return orbit(reflections(W), Identity(W)^W, onInvolutionClasses);
end;
