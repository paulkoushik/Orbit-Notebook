#############################################################################
##
#A  tree.g                                                               GAP4
#B    by Götz Pfeiffer <goetz.pfeiffer@nuigalway.ie>
##
#C  A tree, together with BFS and DFS
##

##  a small tree,  visited by BFS and DFS
##
##  1 -- 3 -- 4 -- 5 -- 6 -- 7
##            |
##            2
##
##  with root 6
##

##  the nodes
nodes := List([1..7], i-> rec(id := i, next := []));

##  the edges
parent := [3,4,4,5,6,0,6];
for i in [1..Length(parent)] do
    if parent[i] > 0 then
        Add(nodes[parent[i]].next, nodes[i]);
    fi;
od;

##  BFS
BFS := function(x, visit)
    local   Q,  y,  z;
    Q:= [x];
    for y in Q do
        visit(y);
        for z in y.next do
            Add(Q, z);
        od;
    od;
end;

##  DFS
DFS := function(x, visit)
    local   z;
    visit(x);
    for z in x.next do
        DFS(z, visit);
    od;
end;

##  a visitor
print := function(x)  Print(x.id, ", ");  end;

##  print as tree
t_indent := 0;  t_nl := false;
t_print := function(x)
    local   c;
    if t_nl then
        Print(RepeatedString(" ", t_indent));  t_nl := false;
    fi;
    Print("-", x.id);
    if Length(x.next) = 0 then
        Print("\n");  t_nl := true;
    fi;
    t_indent := t_indent + 2;
    for c in x.next do  t_print(c);  od;
    t_indent := t_indent - 2;
end;
