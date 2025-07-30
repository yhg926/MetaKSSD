use diagnostics;
use warnings;
if(@ARGV!=2){
 die "*.pl <nodes.dmp> <names.dmp>";
}
open $namef,$ARGV[1] || die "can't open $ARGV[1]:$!";

while(<$namef>){
	chomp;
	next if !/scientific name/;
	($tid,$name)=(split /\t+/)[0,2];
	if(exists $hash{$tid}) {

		$hash{$tid} = $name if length($hash{$tid}) > length($name) ;
	}
	else{
		$hash{$tid} = $name;
	}
}
close $name;

open $node,$ARGV[0] || die "can't open $ARGV[0]:$!";
	while(<$node>){
		chomp;
		($tid,$pid,$rank)=(split /\t+/)[0,2,4];

		print $tid,"\t", $rank,"\t",$pid,"\t",$hash{$tid},"\n";

	}
close $node;

