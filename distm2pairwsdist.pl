use warnings;
use diagnostics;
if(@ARGV!=1){
	die "*.pl <nxn distance matrix> >pairwise distance";
	exit(0);
}
open $dm,$ARGV[0]||die "can't open $ARGV[0]:$!";
$count=0;
while(<$dm>){
	next if /^[^A-Za-z\s]+\n/;
	$sp[$count]=(split /\s/)[0];
	$count++;
}
print "total $count genomes\n";
seek $dm, 0, 0;
while(<$dm>){
	next if /^[^A-Za-z\s]+\n/;
	@dist=(split /\s/);
	$spn = shift @dist;
	if(@dist!=$count){
		print "sp count not match:",scalar @dist,"\t",$count,"\n";
		exit(0);
	}
	for($i=0;$i<$count;$i++){
		print $spn,"\t",$sp[$i],"\t",$dist[$i],"\n";
	}
}
