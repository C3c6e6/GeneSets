#!/usr/bin/env perl

use strict;
use warnings;

#===============================================================================
# Subroutines
#===============================================================================
sub read_file_data {
    my $input_file = "files.txt";
    my $data_dir = "peak_files";
    my %file_data;
    open(IN, $input_file) || die "Could not open $input_file: $!\n";
    while (<IN>) {
        chomp;
        my %record;
        my ($filename, $metadata) = split(/\t/);
        foreach my $pair_string (split(/;\s*/, $metadata)) {
            $pair_string =~ m/(.+)\s*=\s*(.+)/;
            $record{$1} = $2;
        }
        $file_data{"$data_dir/$filename"} = \%record;
    }
    close(IN);
    return(%file_data);
}

sub read_gene_starts {
    my %gene_starts;
    my $max_distance = $_[0];
    my @annotation_files = glob("*.gp.gz");
    my $input_file = $annotation_files[0];
    (scalar(@annotation_files) > 1) && die "Multiple annotation files found!\n";
    open(IN, "gunzip -c $input_file |") || die "Could not open $input_file: $!\n";
    while(<IN>) {
        chomp;
        my ($id, $chrom, $strand, $start, $end) = split(/\t/);
        my $position = $strand eq "+" ? $start : $end;
        my $block = int($position/$max_distance);
        $gene_starts{$chrom}{$block}{$strand}{$id}= $position;
    }
    close(IN);
    return(%gene_starts);
}

sub map_peaks {
    my ($file_data, $gene_starts, $max_distance, %antibody_data) = @_;
    my $c = 0;
    my %targets;
    print STDERR "Mapping peaks ";
    for my $file (keys %$file_data) {
        my $antibody = $$file_data{$file}{"antibody"};
        defined($antibody) || next;
        my $tf = $antibody_data{$antibody}{"target"};
        (-e $file) || next; #Some files listed in files.txt file are missing.
        ($tf =~ m/^eGFP-/) && next; #Ignore e-GFP tagged transcription factors.
        open(IN, "gunzip -c $file |") || die "Could not read $file: $!\n";
        while(<IN>) {
            chomp;
            my @data = split(/\t/);
            my ($chrom, $start, $q, $peak) = @data[0, 1,8,9];
            ($q < 2) && next;
            my $id = find_downstream_transcript_start($chrom, $start+$peak, 
                    $max_distance, $gene_starts);
            defined($id) && push(@{$targets{$tf}}, $id);
        }
        close IN;
        $c++;
        ($c % 10 == 0) && print STDERR ".";
    }
    print STDERR "\n";
    return(%targets);
}

sub find_downstream_transcript_start {
    my ($chrom, $position, $max_distance, $gene_starts) = @_;
    my $closest_distance = $max_distance/2;
    my $closest_id;
    my $block = int($position/$max_distance);
    my %block_gene_starts;
    foreach my $b ($block-1 .. $block+1) {
        exists($$gene_starts{$chrom}{$b}) || next;
        my %this_block = %{$$gene_starts{$chrom}{$b}};
        @block_gene_starts{keys %this_block} = values %this_block;
    }
    for my $strand (("+", "-")) {
        my $sign = "${strand}1";
        for my $id (keys %{$block_gene_starts{$strand}}) {
            my $start = $block_gene_starts{$strand}{$id};
            my $distance = ($start - $position)*$sign;
            if (abs($distance) < $closest_distance) {
                $closest_distance = $distance;
                $closest_id = $id;
            }
        }
    }
    return($closest_id);
}

sub read_conversion_table {
    my $input_file = "conversion_table.txt";
    my %conversion;
    open(IN, $input_file) || die "Could not open $input_file: $!\n";
    <IN>; #Ignore the first line
    while (<IN>) {
        chomp;
        my ($input,$output)=split(/\t/);
        $conversion{$input}=$output;
    }
    close IN;
    return %conversion;
}

sub read_antibody_data {
    my ($input_file, $record);
    my %antibody_records;
    $input_file = "cv.ra";
    $record = {};
    open(IN, $input_file) || die "Could not open $input_file: $!\n";
    while (<IN>) {
        m/^#/ && next;
        chomp;
        s/\s+$//;
        if (length($_) == 0) {
            if (exists($$record{"type"}) && $$record{"type"} eq "Antibody") {
                $antibody_records{$$record{"term"}} = $record;
            }
            $record = {};
            next;
        }
        m/(\w+)\s+(.+)/;
        $$record{$1} = $2;
    }
    close IN;
    return (%antibody_records);
}

sub get_tf_descriptions {
    my %antibody_data = @_;
    my %descriptions;
    foreach my $record (values %antibody_data) {
        (exists($$record{"target"}) && exists($$record{"targetDescription"}))
                || next;
        $descriptions{$$record{"target"}} = $$record{"targetDescription"};
    }
    return(%descriptions);
}

#===============================================================================
# Main
#===============================================================================
my $max_distance = $ARGV[0];

my %file_data = &read_file_data();
my %antibody_data = &read_antibody_data();
my %tf_description = &get_tf_descriptions();
my %gene_starts = &read_gene_starts($max_distance);
my %targets = &map_peaks(\%file_data, \%gene_starts, $max_distance,
        %antibody_data);
my %conversion = &read_conversion_table();

print STDERR "Writing output...";
print STDOUT join("\t", qw(geneID termID termName dbName description))."\n";
foreach my $tf (keys %targets) {
    my $term_id = "ENCODE:$tf";
    my $name = "Targets of transcription factor $tf";
    my $description = "All genes that are shown by ChIP-seq data to have ".
            "transcription factor $tf binding within ".int($max_distance/2).
            " of its transcription start site.";
    exists($tf_description{$tf}) &&
            ($description.="$tf is described as \"$tf_description{$tf}\".");
    my %genes_seen;
    foreach my $target_id (@{$targets{$tf}}) {
        $target_id =~ s/\.\d+$//;
        exists($conversion{$target_id}) || next;
        my $gene_id = $conversion{$target_id};
        $genes_seen{$gene_id} && next;
        $genes_seen{$gene_id} = -1;
        print STDOUT join("\t", $gene_id, $term_id, $name, "ENCODE",
                $description)."\n";
    }
}

print STDERR "\n";
