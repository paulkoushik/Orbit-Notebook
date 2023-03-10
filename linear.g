spinning :=  function(aaa, x, under)
    local list,  a,  y,  z;
    list := [x];
    for y in list do
        for a in aaa do
            z := under(y, a);
            if not z in VectorSpace(Rationals, list) then
                Add(list, z);
            fi;
        od;
    od;
    return list;
end;

sparseVec:= function(vec)
    local   poss;
    poss:= PositionsProperty(vec, v -> v <> Zero(v));
    return rec(pos := poss,  val:= vec{poss});
end;

denseVec:= function(l, v)
    local   vec;
    vec := 0*[1..l];  vec{v.pos} := v.val;
    return vec;
end;

spinning_with_images :=  function(aaa, x, under)
    local   list,  images,  i,  y,  k,  a,  z,  v;
    list := [x];  images := List(aaa, x-> []);  i := 0;
    while i < Length(list) do
        i := i+1;  y := list[i];
        for k in [1..Length(aaa)] do
            a := aaa[k];  z := under(y, a);
            v := SolutionMat(list, z);
            if v = fail then
                Add(list, z);
                v := rec(pos := [Length(list)], val := [1]);
            else
                v := sparseVec(v);
            fi;
            images[k][i]:= v;
        od;
    od;
    return rec(list := list, images := images);
end;

wordTab := function(tab)
    local   word,  i;
    word := [];
    for i in [1..Length(tab)] do
        word{tab[i]} := 0*tab[i] + i;
    od;
    return word;
end;

orbit_with_schreier:= function(aaa, x, under)
    local   list,  reps,  i,  images,  y,  k,  a,  z,  l;
    list := [x];  reps := [()];  i := 0;
    images := List(aaa, x -> []);
    while i < Length(list) do
        i := i+1;  y := list[i];
        for k in [1..Length(aaa)] do
            a := aaa[k];  z := under(y, a);
            l := Position(list, z);
            if l = fail then
                Add(list, z);
                Add(reps, reps[i] * a);
                l := Length(list);
            fi;
            images[k][i] := rec(pos := [l], val := [reps[i] * a / reps[l]]);
        od;
    od;
    return rec(list := list, images := images);
end;

