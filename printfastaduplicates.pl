#print fasta file duplicates

use Bio::SeqIO;
use List::Util 'shuffle';
$in=Bio::SeqIO->new(-file=> '/home/tarakaramji/isomiRs/Testingdatasets/miRNAs_mature_human_miRbaseV21_highconfidence_SNPs.fa', '-format' => 'Fasta');
while(my $seq = $in-> next_seq()){
  $id="";$id=$seq->id;
  $id=~s/\,//g;
  $sequence=$seq->seq();
  $sequence=~tr/[a-z]/[A-Z]/;
for($i=1;$i<=100;$i++)
{
 #     $shuffle="";
  #    @bs=();@bs=split(//,$hash{$mi}{$id});
   #   $shuffle=join("",shuffle(@bs));
      if($sequence)
    {
	#$shufseq="";
	$shufseq=$sequence;
	$shuffle=~tr/[a-z]/[A-Z]/;
	$shufseq=~s/$hash{$mi}{$id}/$shuffle/;
	print ">$i$id\n$shufseq\n";
#	print "$mi\t$id\t$hash{$mi}{$id}\t$shuffle\n";
   }else{
#	print "$mi\t$id\t$hash{$mi}{$id}\n";
      }
    }
  }
# $seq1="";$seq1=substr ($sequence,1,200);

