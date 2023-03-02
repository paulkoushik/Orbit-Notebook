#############################################################################
##
#A  sims.g                                                               GAP4
#B    by Götz Pfeiffer <goetz.pfeiffer@nuigalway.ie>
##
#C  A simple minded implementation of the Schreier-Sims algorithm
##

DeclareAttribute("Sims", IsPermGroup);

is_element_sims := function(group, a)
    local   sims,  pos;
    if IsTrivial(group) then  return a = ();  fi;
    sims := Sims(group);
    pos := Position(sims.list, sims.list[1]^a);
    if pos = fail then  return false;  fi;
    return is_element_sims(sims.stab, a/sims.reps[pos]);
end;

group_plus_one := function(group, a)
    if is_element_sims(group, a) then  return group;  fi;
    return GroupWithGenerators(Concatenation(GeneratorsOfGroup(group), [a]));
end;

orbit_sims := function(aaa, x, under)
    local   list,  reps,  stab,  i,  k,  l,  z;
    list := [x];  reps := [()];  stab := Group([], ());  i := 0;
    while i < Length(list) do
        i := i+1;
        for k in [1..Length(aaa)] do
            z := under(list[i], aaa[k]);  l := Position(list, z);
            if l = fail then
                Add(list, z);
                Add(reps, reps[i] * aaa[k]);
            else   # x^(reps[i] * a) = x^reps[l]
                stab := group_plus_one(stab, reps[i] * aaa[k] / reps[l]);
            fi;
        od;
    od;
    return rec(list := list, reps := reps, stab := stab);
end;

InstallMethod(Sims, "for perm groups", true, [IsPermGroup], 0, function(group)
    return orbit_sims(GeneratorsOfGroup(group), LargestMovedPoint(group), OnPoints);
end);

size_sims := function(group)
    if IsTrivial(group) then  return 1;  fi;
    return Length(Sims(group).list) * size_sims(Sims(group).stab);
end;

random_sims := function(group)
    if IsTrivial(group) then  return ();  fi;
    return random_sims(Sims(group).stab) * Random(Sims(group).reps);
end;
