#############################################################################
##
#A  syt.g                                                                GAP4
#B    by Götz Pfeiffer <goetz.pfeiffer@nuigalway.ie>
##
#C  Partitions and Standard Young Tableaus as orbits under suitable actions
##

#############################################################################
##
##  newtonSum(list) and newtonDif(list)
##
##  Newton's partial sums and forward differences of a list
##
##  properties:
##    newtonDif(newtonSum(list)) = list
##    newtonSum(newtonDif(list)) = list
##

## [l_1, l_2 ..  l_n] -> [l_1, l_2 - l_1, l_3 - l_2 .. l_n -l_{n-1}]
newtonDif := function(l)
    if l = [] then  return [];  fi;
    return Concatenation([l[1]], newtonDif(l{[2..Length(l)]}-l[1]));
end;

## [s_1, s_2 ..  s_n] -> [s_1, s_1 + s_2, s_1 + s_2 + s_3 .. s_1 + ... + s_n]
newtonSum := function(l)
    if l = [] then  return [];  fi;
    return Concatenation([l[1]], newtonSum(l{[2..Length(l)]})+l[1]);
end;

##  right handed versions:
##    newtonDifR(newtonSum(list)) = list
##    newtonSumR(newtonDif(list)) = list
##

## [l_1, l_2 ..  l_n] -> [l_2 - l_1 .. l_n - l_{n-1}, l_n]
newtonDifR := function(l)
    local  ll;
    if l = [] then  return [];  fi;
    ll := Length(l);
    return Concatenation(newtonDifR(l{[1..ll-1]}-l[ll]), [l[ll]]);
end;

## [s_1, s_2 ..  s_n] -> [s_1 + ... + s_n, s_2 + ... + s_n .. s_n]
newtonSumR := function(l)
    local  ll;
    if l = [] then  return [];  fi;
    ll := Length(l);
    return Concatenation(newtonSumR(l{[1..ll-1]})+l[ll], [l[ll]]);
end;

#############################################################################
##
##  Compositions (of n) vs Subsets (of {1..n-1})
##
##  how to turn a subset of [1..n-1] into a composition of n:
##  * find the complement cmp of L in [1..n]
##  * compute Newton differences of cmp
##  e.g. n = 9, L =    45 7  \subseteq [1..8],
##            cmp = 123  6 89,
##            com = 111  3 21.
##
compositionSubset := function(n, set)
    return newtonDif(Difference([1..n], set));
end;

##  how to turn a composition of n into a subset of [1..n-1]
subsetComposition := function(n, com)
    return Difference([1..n], newtonSum(com));
end;

##  examples of orbits:

##  the usual orbit algorithm, again.
orbit := function(aaa, x, under)
    local a, y, z, list;
    list := [x];
    for y in list do
        for a in aaa do
            z := under(y, a);
            if z <> y and not z in list then Add(list, z); fi;
        od;
    od;
    return list;
end;

##  take-away action
takeAway := function(set, s)
    return Difference(set, [s]);
end;

##  power set as orbit under take-away
subsets := function(set)
    return orbit(set, set, takeAway);
end;

##  partitions of n as takeAway orbit of composition classes
partitions := function(n)
    local   takeAwayPartition;
    takeAwayPartition := function(com, a)
        com := compositionSubset(n, takeAway(subsetComposition(n, com), a));
        SortBy(com, x -> -x);
        return com;
    end;
    return Reversed(orbit([1..n-1], [n], takeAwayPartition));
end;

#############################################################################
##
##  SYT: a standard Young tableau is a shortest path in the Young lattice
##

##  normal case action
remove1Hook2 := function(lambda, k)
    local   new;
    new := ShallowCopy(lambda);
    new[k] := new[k] - 1;
    return new;
end;

##  special case action
remove1Hook1 := function(lambda, k)
    if lambda[k] = 1 then
        return lambda{[1..k-1]};
    fi;
    return remove1Hook2(lambda, k);
end;

##  Young lattice and shortest paths.  2-step process.
##  1. build graph "bottom up" from lambda; for each edge, record the row/col
##     positions it corresponds to.
##  2. find shortest paths "top down" from 0
##  (consider setting this up as starting with a list of partitions of n)
##
youngLattice := function(lambda)
    local   grow,  list,  next,  i,  x,  l,  dif,  k;

    ## how to add to the list
    grow := function(list, y, r, c)
        local   pos;
        pos := Position(list, y);
        if pos = fail then
            Add(list, y);  Add(next, []);
            pos := Length(list);
        fi;
        Add(next[i], rec(pos := pos, row := r, col := c));
    end;

    ## orbit-with-edges
    list := [lambda];  next := List(list, x -> []);
    i := 0;
    while i < Length(list) do
        i := i + 1;  x := list[i];  l := Length(x);
        if l > 0 then
            dif := newtonDifR(x);  # ;-)
            grow(list, remove1Hook1(x, l), l, x[l]);  # dif[l] > 0
            for k in Reversed([1..l-1]) do
                if dif[k] > 0 then
                    grow(list, remove1Hook2(x, k), k, x[k]);
                fi;
            od;
        fi;
    od;
    return rec(list := list, next := next);
end;

##  DFS with a visitor that makes shortest paths if unknown.
pathfinder := function(x, next, path)
    local   z,  p;
    if not IsBound(path[x]) then
        path[x] := [];
        for z in next[x] do
            for p in pathfinder(z.pos, next, path) do
                Add(path[x], Concatenation(p, [[z.row, z.col]]));
            od;
        od;
        if path[x] = []  then  Add(path[x], []);  fi;  # empty path
    fi;
    return path[x];
end;

## list of SYTs of shape lambda
standardYTs := function(lambda)
    return pathfinder(1, youngLattice(lambda).next, []);
end;

##  tableau from path
tableau_path := function(path)
    local   tab,  i,  r,  c;
    tab := [];
    for i in [1..Length(path)] do
        r := path[i][1];  c := path[i][2];
        if c = 1 then  tab[r] := [];  fi;
        tab[r][c] := i;
    od;
    return tab;
end;
