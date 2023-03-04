#############################################################################
##
#A  specht.g                                                             GAP4
#B    by Götz Pfeiffer <goetz.pfeiffer@nuigalway.ie>
##
#C  A Specht module example
##

Read("coxeter.g");
Read("linear.g");

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

# how to convert a vector into a list of (pos, val)-pairs
sparseVec:= function(vec)
    local   poss;
    poss:= PositionsProperty(vec, v -> v <> Zero(v));
    return rec(pos := poss,  val:= vec{poss});
end;

# how to convert (pos, val) into a vector of length l
denseVec:= function(l, v)
    local   vec;
    vec := 0*[1..l];  vec{v.pos} := v.val;
    return vec;
end;


# linear orbit with matrices (cf. orbit_with_images)
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

mat_list:= function(list)
    return List(list, v -> denseVec(Length(list), v));
end;

## representing matrices
vvv:= spinning_with_images(perms, vec, Permuted);
mats := List(vvv.images, mat_list);

## Exercise: compute the character ...
