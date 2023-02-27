#############################################################################
##
##  Variants of Relations  (Lists Version.)
##
##
VariantsRelations := function(genrel)
    local   inv,  abs,  variants,  list,  relator,  l,  i,  s,  word,  t;

    inv := word -> -Reversed(word);

    abs := function(s)
        if s < 0 then  return genrel.invr[-s];  else  return s;  fi;
    end;

    variants := List(genrel.gens, x -> []);
    for list in genrel.rels do
        relator :=  Concatenation(list[1], inv(list[2]));
        l := Length(relator);
        for i in [1..l] do
            s := relator[i];    #  u s v = 1, s = (v u)^{-1},  s^-1 = v u
            word := Concatenation(relator{[i+1..l]}, relator{[1..i-1]});
            if s < 0 then
                s := -s;
            else
                word := inv(word);
            fi;
            AddSet(variants[s], List(word, abs));
            t := genrel.invr[s];
            AddSet(variants[t], List(inv(word), abs));
        od;
    od;

    # sort by length.
    for list in variants do
        Sort(list, function(a, b) return Length(a) < Length(b); end);
    od;

    return variants;
end;
