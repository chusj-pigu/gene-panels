## Gene Panels
Repository to store gene panels used for adaptive sampling sequencing projects.

To remove xxx_chr0X_ in front of gene names (used in order to facilitate ordering of genes for coverage plots outputed by MinKnow):

```
awk 'BEGIN {OFS="\t"} {gsub(/[0-9]+_chr[0-9XY]+_/,"",$4); print}'
```