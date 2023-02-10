Read("orbits.g");

isGroupElement:= function(group, a)
    local   x,  orb,  pos,  stab;
    if group = TrivialSubgroup(group) then return a = Identity(group); fi;
    x:= LargestMovedPoint(group);
    orb:= orbit_with_stabilizer(GeneratorsOfGroup(group), x, OnPoints);
    pos:= Position(orb.list, OnPoints(x, a));
    Print(pos, " \c");
    if pos = fail then  return false;  fi;
    stab:= Subgroup(group, Difference(orb.stab, [()]));
    Print("[", Length(GeneratorsOfGroup(stab)), "]\c");
    return isGroupElement(stab, a/orb.reps[pos]);
end;