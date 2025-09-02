#!/usr/bin/env perl
use strict;
use warnings;
use diagnostics;

# Usage: script.pl <kssd2out> <data/gtdbr*_psid2krona_taxonomy.tsv[.gz]>
die "*.pl <MetaKSSD_out> <data/gtdbr*_psid2krona_taxonomy.tsv[.gz]>\n" if @ARGV != 2;

my @ranks = qw(d p c o f g s);
my ($kssd_out_path, $map_path) = @ARGV;

# Open plain or .gz transparently
sub open_maybe_gz {
    my ($path) = @_;
    my $fh;
    if ($path =~ /\.gz$/i) {
        # Use system gzip (no extra Perl modules required)
        open $fh, "-|", "gzip", "-dc", "--", $path
          or die "can't open gzip stream for $path: $!";
    } else {
        open $fh, "<", $path
          or die "can't open $path: $!";
    }
    return $fh;
}

# --- Load psid -> lineage map ---
my %lineage_of;            # psid => "Kingdom|Phylum|...|Species"
my $map_fh = open_maybe_gz($map_path);
while (my $line = <$map_fh>) {
    chomp $line;
    next if $line eq '';
    my ($psid, @rest) = split /\t/, $line, -1;
    $lineage_of{$psid} = join('|', @rest);
}
close $map_fh;

# --- Parse kssd2out and accumulate ---
my %hash;                  # $hash{$smp}[$i]{$rank__taxon} = abundance
my %path_of_taxon;         # rank__taxon => cumulative path "d__..|p__..|..."
my $out_fh = open_maybe_gz($kssd_out_path);

while (my $line = <$out_fh>) {
    chomp $line;
    next if $line eq '';
    my ($smp, $name, $ab) = split /\t/, $line, -1;

    # psid is prefix before first underscore in $name
    my ($psid) = split /_/, $name, 2;

    # sample ID = last path component of $smp
    ($smp) = $smp =~ m{([^/]+)$};

    unless (defined $lineage_of{$psid}) {
        print "$psid is not defined!\n";
        exit(1);
    }

    my @tmp_lineage = split /\|/, $lineage_of{$psid}, -1;

    my $path = '';
    for (my $i = 0; $i < @tmp_lineage && $i < @ranks; $i++) {
        my $rank__taxon = $ranks[$i] . "__" . $tmp_lineage[$i];
        $path = ($path eq '') ? $rank__taxon : ($path . '|' . $rank__taxon);
        $path_of_taxon{$rank__taxon} = $path;
        $hash{$smp}[$i]{$rank__taxon} += $ab;
    }
}
close $out_fh;

# --- Output ---
print "SampleID\tTaxonomy\tRelative_abundance\n";
for my $smp (keys %hash) {
    for my $i (0 .. 6) { # ranks d..s
        next unless defined $hash{$smp}[$i];
        my @sortedtaxa = sort {
            $hash{$smp}[$i]{$b} <=> $hash{$smp}[$i]{$a}
        } keys %{ $hash{$smp}[$i] };

        for my $taxon (@sortedtaxa) {
            my $path = $path_of_taxon{$taxon} // $taxon; # fallback
            print $smp, "\t", $path, "\t", $hash{$smp}[$i]{$taxon}, "\n";
        }
    }
}

