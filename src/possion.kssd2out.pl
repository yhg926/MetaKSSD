use warnings;
use diagnostics;
if(@ARGV != 2){
	die "*.pl <kssdcomposite_output.tsv> <shared_K>";
}
$shkm_thr = defined $ARGV[1] ? $ARGV[1] : 6 ; 

@poissonN=(0,0.02,0.08,0.69,1.57,2.49,3.41,4.31,5.20,6.08,6.94);
$poissonN_thr = 5;
$avgpct9899_offset = 3;

#use origin kssd composite output
#format:Qry\tRef\tShare_kmc\tAvg_kmc\tAvgpct9899_kmc\tMedian_kmc\tTop_kmc
open $kssdf, $ARGV[0] || die "cant open $ARGV[0]: $!";

while(<$kssdf>){
	chomp;
	($sample,$ref, $shkm, $avgpct9899)= (split /\t+/)[0,1,2,4];
	$sample=~s/[^0-9a-zA-Z_.]/_/g;
	#$psid=(split /_/, $ref)[0];
	if( $shkm > $shkm_thr){				
		$depth = $avgpct9899 > $poissonN_thr ? $avgpct9899 - $avgpct9899_offset : $poissonN[int($avgpct9899)];
		
		$data{$sample}->{$ref} = $depth;
		$sum{$sample} += $depth;
	}
}

close $kssdf;

foreach $sample (sort keys %data){
	  foreach $ref ( sort{ $data{$sample}->{$b} <=> $data{$sample}->{$a} } keys %{$data{$sample}}) {
			print $sample,"\t",$ref,"\t",	$data{$sample}->{$ref}/$sum{$sample},"\n";
		}

}
