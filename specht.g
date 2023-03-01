#############################################################################
##
#A  specht.g                                                             GAP4
#B    by Götz Pfeiffer <goetz.pfeiffer@nuigalway.ie>
##
#C  A Specht module example
##

Read("coxeter.g");

N := 5;
gens := transpositions(5);

# a random SYT
tab := [
  [1,4],
  [2,5],
  [3],
];

# how to turn a tableau into a (row) word
wordTab := function(tab)
    local   word,  i,  j;
    word := [];
    for i in [1..Length(tab)] do
        for j in [1..Length(tab[i])] do
            word[tab[i][j]] := i;
        od;
    od;
    return word;
end;

# the (col) tabloid
colWord := wordTab(TransposedMat(tab));
stab := orbit_with_stabilizer(gens, colWord, Permuted).stab;
stab := takeAway(stab, ());
tabs := orbit_with_transversal(stab, tab, OnTuplesTuples);

# action on row words
orb := orbit_with_images(gens, wordTab(tab), Permuted);
words := orb.list;
perms := List(orb.images, PermList);

# the polytabloid of tab
vec := List(words, x-> 0);
for i in [1..Length(tabs.list)] do
    vec[Position(words, wordTab(tabs.list[i]))] := SignPerm(tabs.reps[i]);
od;

# linear orbit
spin :=  function(aaa, x, under)
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

# specht module basis
vvv := spin(perms, vec, Permuted);

## representing matrices
mats := List(perms, p -> List(vvv, v-> SolutionMat(vvv, Permuted(v, p))));

## character ...