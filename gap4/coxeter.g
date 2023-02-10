Read("orbits.g");

cartan_from_graph:= function(n, graph)
    local cartan,  ij,  i,  j;
    cartan:= 2*IdentityMat(n);;
    for ij in graph do
        i:= ij[1];  j:= ij[2];
        cartan[i][j]:= -1;
        cartan[j][i]:= -1;
    od;
    return cartan;
end;

permutation:= function(a, xxx, under)
    return Sortex(List(xxx, x-> under(x, a))) / Sortex(ShallowCopy(xxx));
end;

CoxeterGroup:= function(C)
    local mats,  perms,  phi,  one,  N,  S,  s;
    one:= C^0;  mats:= [];  S:= [1..Length(C)];
    for s in S do
        mats[s]:= C^0;  mats[s]{S}[s]:= one[s] - C[s];
    od;
    phi:= Filtered(orbits(mats, C^0, OnRight), x-> Sum(x) > 0);
    N:= Length(phi);
    Append(phi, -phi);
    perms:= List(mats, m-> Permutation(m, phi, OnRight));
    return rec(perms:= perms, phi:= phi, mats:= mats, N:= N, rank:= Length(S));
end;

graphE8:= [[1,3], [2,4], [3,4], [4,5], [5,6], [6,7], [7,8]];
cartanE8:= 2*IdentityMat(8);;
for ij in graphE8 do
    i:= ij[1]; j:= ij[2]; cartanE8[i][j]:= -1; cartanE8[j][i]:= -1;
od;

CoxeterLength:= function(W, w)
    return Number([1..W.N], i-> i^w > W.N);
end;

permCoxeterWord:= function(W, word)
    if word = [] then return (); fi;
    return Product(W.perms{word});
end;

isLeftDescent:= function(W, w, s)
    return s^w > W.N;
end;

CoxeterWord:= function(W, w)
    local word, a;
    word:= [];
    while w <> () do
        a:= First([1..W.rank], s-> isLeftDescent(W, w, s));
        Add(word, a); w:= W.perms[a] * w;
    od;
    return word;
end;

reducedWord:= function(W, word)
    return CoxeterWord(W, permCoxeterWord(W, word));
end;

longestElement:= function(W, J)
    local   wJ,  s;
    wJ:= ();
    while true do
        s:= First(J, s-> not isLeftDescent(W, wJ, s));
        if s = fail then  return wJ;  fi;
        wJ:= W.perms[s] * wJ;
    od;
end;

prefixes:= function(W, w)
    local onRightDescents;

    onRightDescents:= function(w, s)
        if isLeftDescent(W, w^-1, s) then
            return w * W.perms[s];
        else
            return w;
        fi;
    end;

    return orbit([1..W.rank], w, onRightDescents);
end;

prefixes_with_edges:= function(W, w)
    local onRightDescents;

    onRightDescents:= function(w, s)
        if isLeftDescent(W, w^-1, s) then
            return w * W.perms[s];
        else
            return w;
        fi;
    end;
    
    return orbit_with_edges([1..W.rank], w, onRightDescents);
end;

longestCoset:= function(W, J, L)
    return longestElement(W, J) * longestElement(W, L);
end;

tackOn := function(x, s)
    return Union(x, [s]);
end;

shape := function(W, J)
    local onParabolics;
    onParabolics := function(K, s)
        return OnSets(K, longestCoset(W, K, tackOn(K, s)));
    end;
    return orbit([1..W.rank], J, onParabolics);
end;

takeAway:= function(x, s)
    return Difference(x, [s]);
end;

shapes := function(W)
    local onShapes, S;
    onShapes := function(x, s)
        return Set(shape(W, takeAway(x[1], s)));
    end;
    S := [1..W.rank];
    return orbit(S, shape(W, S), onShapes);
end;

shape_with_edges := function(W, J)
    local onParabolics;
    onParabolics := function(K, s)
        return OnSets(K, longestCoset(W, K, tackOn(K, s)));
    end;
    return orbit_with_edges([1..W.rank], J, onParabolics);
end;


