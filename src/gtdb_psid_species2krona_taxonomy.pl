use warnings;
use diagnostics;
if(@ARGV != 2){
	die "*.pl <psid_GTDBspecies.list> <GTDB_taxnomy.list>";
}

open $db,$ARGV[1] || die "can't open $ARGV[1]:$!";
while(<$db>){
	chomp;
	($spn)=/\;s__([^\;]+)/;
	$hash{$spn} = $_ ;
}
close ($db);

open $spf,$ARGV[0] || die "can't open $ARGV[0]:$!";
while($line= <$spf>){
	chomp $line;

	if($line =~ /^(\d+)_/){ $psid = $1;} 
	else{ die "$line has no psid"};
	$line =~s/^\d+_//g;
	if(!exists $hash{$line}){
		die "$line does not exists in  $ARGV[1]";	
	}

	@ranks = split /;/,$hash{$line} ;
	print $psid;
	for ($i=0;$i< @ranks;$i++){
		$ranks[$i] =~ s/^[dpcofgs]__//;
		print "\t",$ranks[$i];
	}
	print "\n";

}

close $spf;
