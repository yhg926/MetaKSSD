use warnings;
use diagnostics;
if(@ARGV != 2 ){
	die "*.pl <genome.list> <gtdb_genome2tax.tsv>" ;
}

$d = ';'; #rank splitor 
 
open $g2t,$ARGV[1] || die "can't open $ARGV[1]:$!";

$ln = 0;
while(<$g2t>){
	chomp;$ln++;
	($ac,$tax) = (split /\t+/)[0,1];
	if( $ac =~ /(GC[AF]_(\d+))\.\d+/){ #(GC[AF]_\d+\.\d+)
		$ac_num = $2 ;
	}
	
	$sp = (split /$d/, $tax)[6];
	next if !defined $sp;
	if ($sp =~ s/^s__//){
		$sp =~ s/\s/_/g;
	}else{
		print "Warning!: skiping $ARGV[1] line $ln: $ac\t$sp bcs species name do not start with s__\n";
	}		
	$hash{$ac_num} = $sp;

#saw GCF and GCA are equvalent;
#	$ac=~tr/AF/FA/;
#	$hash{$ac} = $sp;	
}
close $g2t;


open $gl,$ARGV[0] || die "can't open $ARGV[0]:$!";

$ln = 0; $taxnum = 0;
while(<$gl>){
	chomp; $ln++;
	@fileds = (split /\t+/);
	$genome = pop @fileds;	

	if( $genome =~ /(GC[AF]_(\d+))\.\d+/){  #(GC[AF]_\d+\.\d+)

		if( exists $hash{$2}){
				$tax = $hash{$2};
				$name2id{$tax} = ++$taxnum if !exists $name2id{$tax} ; 							
				$id = $name2id{$tax};   		
				print $id,"\t",$tax,"\t",$1,"\n";
		}

		else{
			  print 0,"\t","NULL","\t",$1,"\n";
#			print "warning!: $ARGV[0] line $ln: $genome accession was not detected in $ARGV[1]!\n";
#			exit(1);
		}

	}
	else{
		print "aborted!: $ARGV[0] line $ln: $genome contain no accession number!\n";
		exit(1);
	}
}

close $gl;
