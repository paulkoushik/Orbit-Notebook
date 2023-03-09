Read("orbits.g");
Read("variants.g");

orbit_and_more := function(aaa, x, under)
    local   list,  words,  reps,  images,  i,  k,  z,  l;
    list := [x];  words := [[]];  reps := [aaa[1]^0];
    images := List(aaa, x -> []);  i := 0;
    while i < Length(list) do
        i := i+1;
        for k in [1..Length(aaa)] do
            z := under(list[i], aaa[k]);
            l := Position(list, z);
            if l = fail then
                Add(list, z);
                Add(words, onWords(words[i], k));
                Add(reps, reps[i] * aaa[k]);
                l := Length(list);
            fi;
            images[k][i] := l;
        od;
    od;
    return rec(list := list, words := words, reps := reps, images := images);
end;

## Items
ItemFamily := NewFamily("ItemFamily", IsObject);

DeclareRepresentation("IsItem",
    IsComponentObjectRep and IsAttributeStoringRep,
    ["key", "idx", "data", "next"]
);

ItemType := NewType(ItemFamily, IsItem);

Item := function(key)
    local   r;
    r := rec(key := key, data := rec(), next := []);
    return Objectify(ItemType, r);
end;

InstallMethod(PrintObj, "for items", true, [IsItem], 0, function(item)
    Print("Item( ", item!.key, " )");
end);

InstallMethod(\=, "for items", true, [IsItem, IsItem], 0, function(itemL, itemR)
    return itemL!.key = itemR!.key;
end);

InstallMethod(\<, "for items", true, [IsItem, IsItem], 0, function(itemL, itemR)
    return itemL!.key < itemR!.key;
end);

orbit_with_data := function(aaa, item, under)
    local   list,  x,  k,  a,  y,  z;
    list := [item];  item!.idx := 1;  
    item!.data := rec(rep := (), word := []);
    for x in list do
        for k in [1..Length(aaa)] do
            a := aaa[k];
            y := Item(under(x!.key, a));
            z := First(list, z -> z = y);
            if z = fail then
                Add(list, y);  y!.idx := Length(list);
                y!.data := rec(  
                  rep := x!.data.rep * a,
                  word := onWords(x!.data.word, k),
                );
                z := y;
            fi;
            x!.next[k] := z!.idx;
        od;
    od;
    return list;
end;

##  Nodes
NodeFamily := NewFamily("NodeFamily", IsObject);

DeclareRepresentation("IsNode",
  IsComponentObjectRep and IsAttributeStoringRep, 
  ["idx", "word", "data", "next"]
);

NodeType := NewType(NodeFamily, IsNode);

Node := function(word, data)
    local   node;
    node := Objectify(NodeType, rec(word := word, data := data, next := []));
    Add(data.list, node);  node!.idx := Length(data.list);
    data.active := data.active + 1;
    return node;
end;

InstallMethod(PrintObj, "for nodes", true, [IsNode], 0, function(node)
    Print("Node( ", node!.idx, " )");
end);

InstallMethod(String, "for nodes", true, [IsNode], 0, function(node)
    return Concatenation("Node( ", String(node!.idx), " )");
end);

InstallMethod(\=, "for nodes", true, [IsNode, IsNode], 0, function(nodeL, nodeR)
    return nodeL!.idx = nodeR!.idx;
end);

InstallMethod(\<, "for nodes", true, [IsNode, IsNode], 0, function(nodeL, nodeR)
    return nodeL!.idx < nodeR!.idx;
end);

isActive := node -> not IsBound(node!.flat);

flat := function(node)
    while IsBound(node!.flat) do  node := node!.flat;  od;
    return node;
end;

getImage := function(node, s)
    if s < 0 then  s := node!.data.invr[-s];  fi;
    return GetWithDefault(node!.next, s, false);
end;

sprout := function(node, s)
    local   next;
    next := Node(onWords(node!.word, s), node!.data);
    node!.next[s] := next;
    next!.next[node!.data.invr[s]] := node;
    return next;
end;

onNodesPartial := function(node, s)
    local   next;
    next := getImage(node, s);
    if next = false then  return false;  fi;
    return flat(next);
end;

onNodesSprout := function(node, s)
    local   next;
    next := getImage(node, s);
    if next = false then 
        return sprout(node, s);
    else
        return flat(next);  
    fi;
end;

nodeUnderWordSprout := function(node, word)
    local   s;
    for s in word do
        node := onNodesSprout(node, s);
    od;
    return node;
end;

nodeUnderWordPartial := function(node, word)
    local   s;
    for s in word do
        node := onNodesPartial(node, s);
        if node = false then  return node;  fi;
    od;
    return node;
end;

enumerate := function(genrel)
    local  data,  node,  word,  s;

    # initialize.
    data := rec(list := [], active := 0);
    data.invr := genrel.invr;
    data.variants := VariantsRelations(genrel);
    node := Node([], data);
    
    # first close the subgroup tables.
    for word in genrel.sbgp do
        trace(node, word);
    od;
    
    # process nodes in the queue
    for node in data.list do
        for s in genrel.gens do
            process(node, s);
        od;
    od;

    # return data
    return data;
end;

trace := function(node, word)
    local   other;
    other := nodeUnderWordSprout(node, word{[1..Length(word)-1]});
    updateEdge(other, word[Length(word)], node);
end;

process := function(node, s)
    local   variant,  next;
    for variant in node!.data.variants[s] do
        if isActive(node) then
            next := nodeUnderWordPartial(node, variant);
            if next <> false then  updateEdge(node, s, next);  fi;
        fi;
    od;
    if isActive(node) and not IsBound(node!.next[s]) then sprout(node, s); fi;
end;

updateEdge := function(node, s, next)
    setImage(node, s, next);
    setImage(next, node!.data.invr[s], node);
end;

setImage := function(node, s, next)
    local   pair;
    if IsBound(node!.next[s]) then
        pair := Set(List([next, node!.next[s]], flat));
        if Length(pair) = 2 then           # coincidence: stack!
            mergeNodes(pair[2], pair[1]);
        fi;
    else
        node!.next[s] := next;           # deduction!
    fi;
end;

mergeNodes := function(node, other)
    local   s;
    node!.flat := other;
    node!.data.active := node!.data.active - 1;
    for s in PositionsBound(node!.next) do
        updateEdge(other, s, node!.next[s]);
    od;
end;
