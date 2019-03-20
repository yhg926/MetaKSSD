use warnings;
use diagnostics;
die"*.pl <pair_distance.list> <sepcies number>" if @ARGV!=2;
open(IN,$ARGV[0])||die "can't open $ARGV[0]:$!";
#print "enter number of sp \n";
$spn=$ARGV[1];

%hash=();
%namecode=();
$ncode=0;
while($line=<IN>){
chomp $line;
($a,$b,$c)=(split "\t" ,$line)[0,1,2] ;
for $na($a,$b){
    if (! exists $namecode{$na}) {
     $name=sprintf "%03s", ++$ncode ;
     $namecode='Species'.$name;
     $namecode{$na}=$namecode;   
 };

}
$hash{$a}->{$b}=$c ;
$hash{$b}->{$a}=$c ;
#print $a,"\t",$b,"\t",$hash{$a}->{$b},"\n";
};

print $spn;
#print "===",$hash{'E_coli_UTI89'}->{'E_fer_ATCC_35469'},"===" ;
@ar=(keys %hash);
#print scalar @ar,"\n";
foreach $ele (@ar){
	print "\n";
	#printf "%-10s", substr($ele,-10,10); #phlip only accept sample name 9 or 10 in length !
        print $namecode{$ele};
	foreach $ele2(@ar){
#		print $ele2 ;
		print " ";
#		print "***",$hash{$ele}->{$ele2},"^^^^" ;
		if ($ele eq $ele2){ 
			printf "%.10f" , 0 ;
			}
		else{
			printf "%.10f" , exists $hash{$ele}->{$ele2}? $hash{$ele}->{$ele2} : $hash{$ele2}->{$ele} ;
			 ;};
  };	
};
print "\n";

open(OUT, ">./namecode")||die"can't open: $!";
for (keys %namecode){
  print OUT $_ ,"\t", $namecode{$_} ,"\n";
}



