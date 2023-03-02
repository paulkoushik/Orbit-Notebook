# use schreier generators as images
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

swaps:= transpositions(4);
os:= orbit_with_schreier(swaps, 4, OnPoints);
mats:= List(os.images, mat_list);

# need to be able to add and multiply 0 and perms
InstallMethod(\*, "for 0 and perm", true, [IsInt, IsPerm], 0, function(zero, perm)
    return zero;
end);
InstallMethod(\*, "for perm and 0", true, [IsPerm, IsInt], 0, function(perm, zero)
    return zero;
end);
InstallTrueMethod( IsAdditiveElementWithInverse, IsPerm );
InstallMethod(\+, "for 0 and perm", true, [IsInt, IsPerm], 0, function(zero, perm)
    return perm;
end);
InstallMethod(\+, "for perm and 0", true, [IsPerm, IsInt], 0, function(perm, zero)
    return perm;
end);
