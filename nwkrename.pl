use diagnostics;
use warnings;
die"*.pl <outtree/outfile>" if @ARGV != 1;
%code2name=();
open(IN,"<./namecode")||die"can't open ./namecode:$!";
while($line=<IN>){
    chomp $line;
    ($name,$code)=(split /\t/, $line)[0,1] ; 
    $code2name{$code}=$name ;
};

open($nwk,"$ARGV[0]")||die"can't open $ARGV[0] :$!";
$/=undef ;
$line=<$nwk>;
$line=~s/(Species\d\d\d)/$code2name{$1}/g ;
print $line;














