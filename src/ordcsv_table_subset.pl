use warnings;
use diagnostics;
if(@ARGV!=2){
	die "*.pl <all.csv> <selected_1stcolname>";
}

open $csv, $ARGV[0]||die "cant open $ARGV[0]:$!";
$spliter = ',';
while(<$csv>){
    chomp;
    $id = (split /$spliter/)[0];
    $hash{$id} = $_ ;
}
close($csv);


open $slct,$ARGV[1]||die "cant open $ARGV[1]:$!";

while(<$slct>){
	chomp;
	$id=(split /\t/)[0];
	if(!exists $hash{$id}) {
		print "==no match==\n";
	}
	else{
		print "$hash{$id}\n";
	}

}
close($slct);
