use warnings;
use diagnostics;
if(@ARGV!=2){
	die "*.pl <kssd2out> <gtdbr214_psid2krona_taxonomy.tsv>";
}

@ranks=('d','p','c','o','f','g','s');
open $map,$ARGV[1] || die "can't open $ARGV[1]:$!";
while(<$map>){
	chomp;
	@line = (split /\t/);
	$psid = shift @line;
	$lineage[$psid] = join '|', @line;  
}
close $map;

open $out,$ARGV[0] || die "can't open $ARGV[0]:$!"; 
while(<$out>){
	chomp;
	($smp,$name,$ab) = (split /\t/);
	$psid = (split /_/, $name)[0];	

	@tmp = (split /\//,$smp);  
	$smp = pop @tmp ;	

	if( !defined $lineage[$psid] ){
		print "$psid is not defined!\n";
		exit(1);
	}
	@tmp_lineage  = split /\|/, $lineage[$psid];
	$path="";
	for($i=0;$i<@tmp_lineage;$i++){
		$rank__taxon = $ranks[$i]."__".$tmp_lineage[$i]; 
		$path eq ""? $path.= $rank__taxon : $path = $path.'|'.$rank__taxon;	
		$pathash{$rank__taxon} = $path;
		$hash{$smp}->[$i]->{$rank__taxon} += $ab; 
		
	}

}
close $out;

print "SampleID\tTaxonomy\tRelative_abundance\n";
foreach $smp (keys %hash){
	for($i=0; $i <7; $i++){
		@sortedtaxa = sort { $hash{$smp}->[$i]->{$b} <=>  $hash{$smp}->[$i]->{$a}  } keys %{$hash{$smp}->[$i]};
		foreach $taxon (@sortedtaxa){
			print $smp,"\t", $pathash{$taxon},"\t",$hash{$smp}->[$i]->{$taxon},"\n";

		}	

	}

}
