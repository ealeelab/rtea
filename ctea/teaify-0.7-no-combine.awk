#!/bin/awk -f
BEGIN { OFS="\t" }
{
fa="."
AluY=L1HS=SVA=PolyA=PolyT=0
clv="-"
split($1, a, ";")


clv=0

if ( $3 != "*" && !($1 ~ /^@/) ) {

family=fa=""

   for (i=1;i<=NF;i++) {
		if($i ~ /^XA:Z:/) {
			family=$i
			sub(/^XA:Z:/, "", family)
		}
    }

    split(family, famil, ";")
    for (i=1; i <= length(famil); i++) {
        split(famil[i], fam, ",")
        if (fam[1] == "AluY") AluY++
        else if (fam[1] == "L1HS") L1HS++
        else if (fam[1] == "SVA") SVA++
        else if (fam[1] == "PolyA") PolyA++
    }

    pre = 0;
    if (AluY != 0) {
        fa = "AluY:" AluY
        pre = 1
    }
    if (L1HS != 0) {
        if (pre == 1) fa = fa ","
        fa = fa "L1HS:" L1HS
        pre = 1
    }
    if (SVA != 0) {
        if (pre == 1) fa = fa ","
        fa = fa "SVA:" SVA
        pre = 1
    }
    if (PolyA != 0) {
        if (pre == 1) fa = fa ","
        fa = fa "PolyA:" PolyA
    }
    if (length(fa) == 0) {
        fa = "."
    }

    print a[1], a[2], $3, fa, a[3]
}
}
