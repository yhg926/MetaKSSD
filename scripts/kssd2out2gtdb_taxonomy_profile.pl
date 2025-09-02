#!/usr/bin/env perl
use strict;
use warnings;
use diagnostics;

# Usage: script.pl <kssd2out> <data/gtdbr*_psid2krona_taxonomy.tsv[.gz]>
die "*.pl <MetaKSSD_out> <data/gtdbr*_psid2krona_taxonomy.tsv[.gz]>\n" if @ARGV != 2;

my @ranks = qw(d p c o f g s);
my ($kssd_out_path, $map_path) = @ARGV;

sub open_maybe_gz {
    my ($path) = @_;
    my $fh;
    if ($path =~ /\.gz$/i) {
        open $fh, "-|", "gzip", "-dc", "--", $path
          or die "can't open gzip stream for $path: $!";
    } else {
        open $fh, "<", $path
          or die "can't open $path: $!";
    }
    return $fh;
}

# --- trim helpers ---
sub trim { my ($x)=@_; $x =~ s/^\s+//; $x =~ s/\s+$//; return $x; }
sub trim_all { return map { trim($_) } @_; }

# --- Load psid -> lineage map ---
my %lineage_of;  # psid => "Kingdom|Phylum|...|Species"
my $map_fh = open_maybe_gz($map_path);
while (my $line = <$map_fh>) {
    $line =~ s/\r?\n$//;              # chomp plus CRLF safety
    next if $line eq '' || $line =~ /^\s*#/;  # skip blanks/comments
    my ($psid, @rest) = split /\t/, $line, -1;
    ($psid, @rest) = trim_all($psid, @rest);
    next if $psid eq '';              # guard
    $lineage_of{$psid} = join('|', @rest);
}
close $map_fh;

# --- Parse kssd2out and accumulate ---
my %hash;            # $hash{$smp}[$i]{$rank__taxon} = abundance
my %path_of_taxon;   # rank__taxon => full path
my $out_fh = open_maybe_gz($kssd_out_path);

while (my $line = <$out_fh>) {
    $line =~ s/\r?\n$//;
    next if $line eq '';
    my ($smp, $name, $ab) = split /\t/, $line, -1;
    ($smp, $name) = trim_all($smp, $name);
    $ab = 0 + $ab;  # numeric

    # psid is prefix before first underscore in $name
    my ($psid) = split /_/, $name, 2;
    $psid = trim($psid);

    # sample ID = last path component of $smp
    ($smp) = ($smp =~ m{([^/]+)$});

    unless (exists $lineage_of{$psid}) {
        print "$psid is not defined!\n";
        exit(1);
    }

    my @tmp_lineage = split /\|/, $lineage_of{$psid}, -1;
    my $path = '';
    for (my $i = 0; $i < @tmp_lineage && $i < @ranks; $i++) {
        my $tax = trim($tmp_lineage[$i]);
        my $rank__taxon = $ranks[$i] . "__" . $tax;
        $path = ($path eq '') ? $rank__taxon : ($path . '|' . $rank__taxon);
        $path_of_taxon{$rank__taxon} = $path;
        $hash{$smp}[$i]{$rank__taxon} += $ab;
    }
}
close $out_fh;

# --- Output ---
print "SampleID\tTaxonomy\tRelative_abundance\n";
for my $smp (keys %hash) {
    for my $i (0 .. 6) {
        next unless defined $hash{$smp}[$i];
        my @sortedtaxa = sort {
            $hash{$smp}[$i]{$b} <=> $hash{$smp}[$i]{$a}
        } keys %{ $hash{$smp}[$i] };
        for my $taxon (@sortedtaxa) {
            my $path = $path_of_taxon{$taxon} // $taxon;
            print $smp, "\t", $path, "\t", $hash{$smp}[$i]{$taxon}, "\n";
        }
    }
}

