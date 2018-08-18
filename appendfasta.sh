#append two fast sequences from two different files with header :
#!>RNA1
#AATGACGATGACGATGACAGAT
#>RNA2
#ATAGATGGGCAGTAGAGA
#File2:
#
#>mRNA1
#ATGGAGATGAGAT
#>mRNA2
#AGATGGGGATGA
#Ouput file should be

#>RNA1:mRNA1
#AATGACGATGACGATGACAGATATGGAGATGAGAT
#>RNA2:mRNA2
#ATAGATGGGCAGTAGAGAAGATGGGGATGA


paste -d '' file1 file2 | sed 's/>/:/2' or  paste -d '' file1 <(tr ">" ":" <file2)
