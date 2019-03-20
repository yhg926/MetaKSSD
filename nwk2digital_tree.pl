#!/usr/bin/perl -w
#yhg copyright
use warnings;
use diagnostics;
if(@ARGV!=1){
	die "*.pl <*.nwk>";
	exit(1); 
} 
$xoffset=1;
$yoffset=1;
$max_width=100; #100 '-' maximal in one line 
$line_max=128; #maximal character allowed in one line

$set_branch='-';
$set_vert='!';
$set_shoulder='+';

open $infile,$ARGV[0]||die "can not open:$!";
@lines=<$infile>;
chomp @lines;
$nwk = join '',@lines;

$level=$i=0;
@taxa=();
$y=$yoffset;

##begin to compute Y
while($nwk=~/(\()|([^\(\)\:\,\;]+:-?[\.0-9]+)|(\)[\.0-9]*:-?[\.0-9]+)/g){
	if($1){
		$level++;
		$branch[$level]=0; # so hard to debug :$branch{$level}=0
	}
	elsif($2){
		($name,$len)=(split /:/,$2);
		$len=0 if $len < 0;
		push @{$current[$level]},$i ; # current level index 
		$branch[$level]=$len if $branch[$level] < $len;#longest branch at current level
		push @taxa,[$name,"leaf",$len,undef,$y];
		$i++;
		$y+=2;	
	}
	elsif($3){
		($bv,$len)=(split /:/,$3);
		$len=0 if $len < 0;
		$bv=~s/\)//;
		#########
		if (@{$current[$level]} !=2){
			print scalar @{$current[$level]},"error detect\n";
			exit(1);
		}
		my($y1,$y2)=($taxa[$current[$level]->[0]]->[4],$taxa[$current[$level]->[1]]->[4]);
		$nodey=int(($y1+$y2)/2);
		my @sub_nodes = @{$current[$level]};
		if ($bv=~/[0-9.]+/){
			$fbv=sprintf("%02d",$bv);
			$fbv="%%" if $bv >99;
		}else{$fbv=$set_branch.$set_vert};		
		push @taxa, [\@sub_nodes,$fbv,$len,undef,$nodey];	
		@{$current[$level]}=();
		$tmplen=$len + $branch[$level];
####################################		
		$level--;
		$branch[$level]=$tmplen if $branch[$level] < $tmplen;		
		push @{$current[$level]},$i ;		
		$i++;				
	}
}


$longestbr= $branch[$level];
sub fill_x($$);
#fill x coordinates;
foreach (@{$current[$level]}){
	fill_x($_,$xoffset);
}

sub map_coord($$);
foreach (@{$current[$level]}){
        map_coord($_,$xoffset);
}
my($first,$last)=@{$current[$level]}[0,-1];
my($firsty,$lasty)=($taxa[$first]->[4],$taxa[$last]->[4]); 
@{$mapped_tree[$xoffset]}[$firsty..$lasty]=($set_shoulder, ($set_vert) x ($lasty-$firsty-1),$set_shoulder);
#=pod
for(my $m=0;$m<$y+1;$m++){
	for (my $n=0;$n<$line_max;$n++){
		if(!defined $mapped_tree[$n]->[$m]){
			print " ";
		}else{
			print $mapped_tree[$n]->[$m];
		}
	}
	print "\n";	
}

print "\n";
print " " x ($xoffset + 2),'+','-' x 4,'+', "\n";
print " " x ($xoffset + 2), sprintf("%.4f",$longestbr*5/$max_width),"\n";

#=cut

sub fill_x($$){ # 'sub fill_x($$){' will report error related with prototype. 
	my $ind = shift @_;
	my $x = shift @_;
	my $x_base = int($taxa[$ind]->[2]/$longestbr*$max_width) + $x + 1;
	if ($taxa[$ind]->[1] ne "leaf"){ #inner node
		$taxa[$ind]->[3] = $x_base + 1;
		fill_x($taxa[$ind]->[0]->[0],$taxa[$ind]->[3]);
		fill_x($taxa[$ind]->[0]->[1],$taxa[$ind]->[3]);
	}
	else{
		$taxa[$ind]->[3] = $x_base;
	}		
}

sub map_coord($$){
	my $ind=shift @_;
	my $pre_x=shift @_;
	my $x = $taxa[$ind]->[3];	
	for(my $i=$pre_x+1;$i<$x;$i++){
		$mapped_tree[$i]->[$taxa[$ind]->[4]]=$set_branch ;
	}
	if ($taxa[$ind]->[1] ne "leaf"){ #inner node
		my($sub_y1,$sub_y2) = ($taxa[$taxa[$ind]->[0]->[0]]->[4],$taxa[$taxa[$ind]->[0]->[1]]->[4]);
		@{$mapped_tree[$x]}[$sub_y1..$sub_y2]=($set_shoulder, ($set_vert) x ($sub_y2-$sub_y1-1),$set_shoulder);
		
	#	$taxa[$ind]->[1] = "%%" if $taxa[$ind]->[1] > 99;
		$mapped_tree[$x-1]->[$taxa[$ind]->[4]]=substr($taxa[$ind]->[1],0,1);
		$mapped_tree[$x]->[$taxa[$ind]->[4]]=substr($taxa[$ind]->[1],1,1);	
		map_coord($taxa[$ind]->[0]->[0],$x);
		map_coord($taxa[$ind]->[0]->[1],$x);	
	}else{
		$mapped_tree[$x]->[$taxa[$ind]->[4]]=$taxa[$ind]->[0];
	}	
}














