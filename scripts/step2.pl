#!/usr/bin/perl -w
use strict;
use Getopt::Long;

################################################################################
# Script: step2.pl
# Purpose: Generate enhancer-gene candidate pairs according to genomic distance.
#          For each enhancer, genes on the same chromosome with TSS within 1 Mb
#          of the enhancer center are retained as putative target genes.
# Inputs : --enhfile enhancer BED-like file
#          --genefile gene GFF3 annotation file
# Output : --outfile enhancer-gene distance table
################################################################################

# Command-line parameters
my ($enhfile, $genefile, $outfile, $help);
GetOptions(
    "enhfile=s"  => \$enhfile,
    "genefile=s" => \$genefile,
    "outfile=s"  => \$outfile,
    "help!"      => \$help,
);

if (defined $help or not ($enhfile and $genefile and $outfile)) {
    print "Usage: $0 --enhfile <enhancer_bed> --genefile <gene_gff3> --outfile <output_file>\n";
    exit;
}

# 读取 enhancer 文件
my %enhHash;
open(my $ENH, "<", $enhfile) or die "Cannot open $enhfile: $!\n";
while (<$ENH>) {
    chomp;
    my @fields = split /\t/;
    next unless scalar(@fields) >= 4;
    my $enh_id = $fields[0];
    my $chrom  = $fields[1];
    my $start  = $fields[2];
    my $end    = $fields[3];
    my $center = int(($start + $end) / 2);
    $enhHash{$enh_id} = {
        chrom  => $chrom,
        start  => $start,
        end    => $end,
        center => $center
    };
}
close $ENH;

# 读取 gene 文件，提取 gene_id, gene_name, TSS, gene_start, gene_end, strand
my %geneHash;
open(my $GENE, "<", $genefile) or die "Cannot open $genefile: $!\n";
while (<$GENE>) {
    chomp;
    next if /^#/;
    my @fields = split /\t/;
    next unless $fields[2] eq "gene";
    my $chrom  = $fields[0];
    my $start  = $fields[3];
    my $end    = $fields[4];
    my $strand = $fields[6];
    my $attr   = $fields[8];

    my ($geneid, $genename);
    if ($attr =~ /gene_id=([^;]+)/) {
        $geneid = $1;
        $geneid =~ s/\.\d+$//;
    }
    if ($attr =~ /gene_name=([^;]+)/) {
        $genename = $1;
    } else {
        $genename = $geneid;
    }

    my $tss = ($strand eq "+") ? $start : $end;
    $geneHash{$chrom}{$geneid} = {
        name   => $genename,
        tss    => $tss,
        start  => $start,
        end    => $end,
        strand => $strand
    };
}
close $GENE;

# 输出结果
open(my $OUT, ">", $outfile) or die "Cannot write to $outfile: $!\n";
print $OUT "enhancer_id\tenhancer_chr\tenhancer_start\tenhancer_end\tenhancer_center\tgene_id\tgene_name\tTSS\tgene_start\tgene_end\tgene_strand\tdistance\n";

foreach my $enhid (keys %enhHash) {
    my $chrom      = $enhHash{$enhid}->{chrom};
    my $enhstart   = $enhHash{$enhid}->{start};
    my $enhend     = $enhHash{$enhid}->{end};
    my $enhcenter  = $enhHash{$enhid}->{center};
    next unless exists $geneHash{$chrom};

    foreach my $geneid (keys %{$geneHash{$chrom}}) {
        my $genename   = $geneHash{$chrom}{$geneid}->{name};
        my $tss        = $geneHash{$chrom}{$geneid}->{tss};
        my $genestart  = $geneHash{$chrom}{$geneid}->{start};
        my $geneend    = $geneHash{$chrom}{$geneid}->{end};
        my $strand     = $geneHash{$chrom}{$geneid}->{strand};
        my $dist       = abs($enhcenter - $tss);

        if ($dist <= 1_000_000 and $dist > 1000) {
            print $OUT join("\t", 
                $enhid, $chrom, $enhstart, $enhend, $enhcenter,
                $geneid, $genename, $tss, $genestart, $geneend, $strand, $dist
            ), "\n";
        }
    }
}
close $OUT;
