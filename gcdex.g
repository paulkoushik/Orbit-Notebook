#############################################################################
##
#A  gcdex.g                                                              GAP4
#B    by Götz Pfeiffer <goetz.pfeiffer@nuigalway.ie>
##
#C  Extended Euclid's algorithm by DFS.
##  Input: integers a, b.  Output: x, y such that gcd(a, b) = ax + by.
##

##  recursive version.
gcdex0:= function(a, b)
    local   qr,  st;
    if b = 0 then return [1, 0]; fi;
    qr:= QuotientRemainder(a, b);
    st:= gcdex0(b, qr[2]);
    return [st[2], st[1] - qr[1] * st[2]];
end;

## 2 visitor version of DFS
DFS2:= function(x, visit1, visit2)
    local   val,  z;
    visit1(x);
    for z in x.next do
        DFS2(z, visit1, visit2);
        visit2(x);
    od;
end;

##  explicit nodes version.
##  each node eventually has fields ab and xy such that gcd(a, b) = ax + by.
gcdex:= function(a, b)
    local   node,  rise,  fall;
    node:= rec(ab:= [a, b]);

    # extend tree if necessary.
    rise:= function(node)
        node.xy:= [1, 0];  node.next:= [];
        if node.ab[2] <> 0 then
            node.qr:= QuotientRemainder(node.ab[1], node.ab[2]);
            Add(node.next, rec(ab:= [node.ab[2], node.qr[2]]));
        fi;
    end;

    # propagate xy down the tree.
    fall:= function(node)
        local   xy;
        xy:= node.next[1].xy;
        node.xy:= [xy[2], xy[1] - node.qr[1] * xy[2]];
    end;

    DFS2(node, rise, fall);

    return node.xy;
end;
