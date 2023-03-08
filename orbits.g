#############################################################################
##
#A  orbits.g                                                             GAP4
#B    by Götz Pfeiffer <goetz.pfeiffer@nuigalway.ie>
##
#C  Examples of Orbit Algorithms and Applipcations.
##

##  a list of transpositions
transpositions := n -> List([1..n-1], j-> (j,j+1));

##  actions
onPoints := function(x, a)  return x^a;  end;

onRight := function(x, a)  return x * a;  end;

onClasses := function(x, a)
    return OnRight(Representative(x), a)^ActingDomain(x);
end;

onGroups := function(x, a)
    return ClosureGroup(x, a);
end;

onSubgroupClasses := function(x, a)
    return onGroups(Representative(x), a)^ActingDomain(x);
end;

##  orbits
orbit := function(aaa, x, under)
    local a, y, z, list;
    list := [x];
    for y in list do
        for a in aaa do
            z := under(y, a);
            if z <> y and not z in list then  Add(list, z);  fi;
        od;
    od;
    return list;
end;

orbits := function(aaa, xxx, under)
    local a, y, z, list;
    list := ShallowCopy(xxx);
    for y in list do
        for a in aaa do
            z := under(y, a);
            if z <> y and not z in list then  Add(list, z);  fi;
        od;
    od;
    return list;
end;

onWords := function(word, s)
    return Concatenation(word, [s]);
end;

orbit_with_words := function(aaa, x, under)
    local   list,  words,  i,  k,  z;
    list := [x];  words := [[]];  i := 0;
    while i < Length(list) do
        i := i+1;
        for k in [1..Length(aaa)] do
            z := under(list[i], aaa[k]);
            if not z in list then
                Add(list, z);
                Add(words, onWords(words[i], k));
            fi;
        od;
    od;
    return rec(list := list, words := words);
end;

orbit_with_transversal := function(aaa, x, under)
    local   list,  reps,  i,  k,  z;
    list := [x];  reps := [aaa[1]^0];  i := 0;
    while i < Length(list) do
        i := i+1;
        for k in [1..Length(aaa)] do
            z := under(list[i], aaa[k]);
            if not z in list then
                Add(list, z);
                Add(reps, reps[i] * aaa[k]);
            fi;
        od;
    od;
    return rec(list := list, reps := reps);
end;

orbit_with_stabilizer := function(aaa, x, under)
    local   list,  reps,  stab,  i,  k,  l,  z;
    list := [x];  reps := [aaa[1]^0];  stab := [];  i := 0;
    while i < Length(list) do
        i := i+1;
        for k in [1..Length(aaa)] do
            z := under(list[i], aaa[k]);
            l := Position(list, z);
            if l = fail then
                Add(list, z);
                Add(reps, reps[i] * aaa[k]);
            else   # x^(reps[i] * a) = x^reps[l]
                Add(stab, reps[i] * aaa[k] / reps[l]);
            fi;
        od;
    od;
    return rec(list := list, reps := reps, stab := stab);
end;

orbit_with_edges := function(aaa, x, under)
    local   list,  edges,  i,  k,  l,  z;
    list := [x];  edges := [];  i := 0;
    while i < Length(list) do
        i := i+1;
        for k in [1..Length(aaa)] do
            z := under(list[i], aaa[k]);
            l := Position(list, z);
            if l = fail then
                Add(list, z);
                Add(edges, [i, Length(list)]);
            else
                Add(edges, [i, l]);
            fi;
        od;
    od;
    return rec(list := list, edges := edges);
end;

## applications
elements := function(group)
    return orbit(GeneratorsOfGroup(group), Identity(group), OnRight);
end;

class := function(group, x)
    return Set(orbit(GeneratorsOfGroup(group), x, OnPoints));
end;

conjugacyClasses := function(group)
    local gens;
    gens := GeneratorsOfGroup(group);
    return orbit(orbits(gens, gens, OnPoints), Identity(group)^group, onClasses);
end;

subgroups := function(group)
    return orbit(Elements(group), TrivialSubgroup(group), onGroups);
end;

subgroupClasses := function(group)
    return orbit(Zuppos(group), TrivialSubgroup(group)^group, onSubgroupClasses);
end;

sizeOfGroup := function(group)
    local x, orb, stab;
    if group = TrivialSubgroup(group) then return 1; fi;
    x := LargestMovedPoint(group);
    orb := orbit_with_stabilizer(GeneratorsOfGroup(group), x, OnPoints);
    stab := Subgroup(group, Difference(orb.stab, [()]));
    return sizeOfGroup(stab) * Length(orb.reps);
end;

randomGroupElement := function(group)
    local x, orb, stab;
    if group = TrivialSubgroup(group) then return Identity(group); fi;
    x := LargestMovedPoint(group);
    orb := orbit_with_stabilizer(GeneratorsOfGroup(group), x, OnPoints);
    stab := Subgroup(group, Difference(orb.stab, [()]));
    return randomGroupElement(stab) * Random(orb.reps);
end;

orbit_with_edges := function(aaa, x, under)
    local   list,  edges,  i,  k,  l,  z;
    list := [x];  edges := [];  i := 0;
    while i < Length(list) do
        i := i+1;
        for k in [1..Length(aaa)] do
            z := under(list[i], aaa[k]);
            l := Position(list, z);
            if l = fail then
                Add(list, z);
                l := Length(list);
            fi;
            Add(edges, [i, l]);
        od;
    od;
    return rec(list := list, edges := edges);
end;

orbit_with_images := function(aaa, x, under)
    local   list,  images,  i,  k,  l,  z;
    list := [x];  i := 0;
    images := List(aaa, x -> []);
    while i < Length(list) do
        i := i+1;
        for k in [1..Length(aaa)] do
            z := under(list[i], aaa[k]);
            l := Position(list, z);
            if l = fail then
                Add(list, z);
                l := Length(list);
            fi;
            images[k][i] := l;
        od;
    od;
    return rec(list := list, images := images);
end;

list_with_index := list -> List([1..Length(list)], i -> [i, list[i]]);
