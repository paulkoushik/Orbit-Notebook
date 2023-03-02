Read("orbits.g");

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
