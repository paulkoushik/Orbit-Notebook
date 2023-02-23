orbit:= function(aaa, x, under)
    local a, y, z, list;
    list:= [x];
    for y in list do
        for a in aaa do
            z:= under(y, a);
            if z <> y and not z in list then Add(list, z); fi;
        od;
    od;
    return list;
end;

orbits:= function(aaa, xxx, under)
    local a, y, z, list;
    list:= ShallowCopy(xxx);
    for y in list do
        for a in aaa do
            z:= under(y, a);
            if z <> y and not z in list then Add(list, z); fi;
        od;
    od;
    return list;
end;

onInvolutions := function(x, s)
    if x * s = s * x then return OnRight(x, s); fi;
    return OnPoints(x, s);
end;

onInvolutionClasses:= function(x, a)
    local y;
    y:= Representative(x);
    if y^a <> y then return x; fi;
    return ConjugacyClass(x.group, OnRight(y, a));
end;

involutionClasses := function(group)
    local gens;
    gens:= Generators(group);
    return orbit(orbits(gens, gens, OnPoints), ConjugacyClass(group, ()), onInvolutionClasses);
end;
