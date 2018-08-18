##cuting first three strings of word in a file
 cut -c-3 test | while read line; do echo $line; done
###
##to replace the annotation. Takes file Annotation change and replaces annotation in tree file
sed 's| *\([^ ]*\) *\([^ ]*\).*|s/\1/\2/g|' <test | sed -f- ./RAxML_bipartitionsBranchLabels.AUTO >SMXL_Raxml_1_12_nameedited
