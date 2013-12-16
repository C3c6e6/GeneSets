#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

#===============================================================================
# Subroutines
#===============================================================================
sub read_organism_file {
    my $file_name = $_[0];
    open(IN, $file_name) || die "Could not open $file_name: $!\n";
    my @organisms = <IN>;
    chomp(@organisms);
    close IN;
    return(@organisms);
}

sub create_symbol2id {
    my $organism_name = $_[0];
    my $script = dirname($0)."/create_conversion_table.R";
    my %conversion;
    $organism_name =~ s/ /./g;
    for my $symbol_type ("SYMBOL", "ALIAS") {
        open(IN, "$script $organism_name $symbol_type ENTREZID 2> /dev/null |")
                || die "Could not open pipe: $!\n";
        <IN>; #Skip first line with column names.
        while (<IN>) {
            chomp;
            my ($symbol, $entrezid) = split(/\t/);
            #undef values indicate ambiguous symbols
            $conversion{$symbol} = (exists($conversion{$symbol}) &&
                    (!defined($conversion{$symbol}) ||
                    $conversion{$symbol} ne $entrezid)) ? undef : $entrezid;
        }        
        close IN;
    }
    return(%conversion);
}

sub read_tf_data {
    my $org = $_[0];
    my $input_file = "${org}_tf.txt";
    my %data;
    open(IN, $input_file) || die "Could not open $input_file: $!\n";
    my @field_names = split(/\t/, <IN>);
    chomp(@field_names);
    while(<IN>) {
        chomp;
        my %record;
        @record{@field_names} = split(/\t/);
        $data{$record{"gene_symbol"}} = \%record;
    }
    close IN;
    return %data;
}

sub process_target_file {
    my ($long_name, $short_name, $package_dir, $input_file, $tf_data,
        $symbol2entrezid, $output_file);
    my %bad_symbols;    
    ($long_name, $short_name, $package_dir, $tf_data, $symbol2entrezid) = @_;
    $input_file = "${short_name}_tf_regulated_gene.txt";
    $long_name =~ s/\s/./g;
    $output_file = "$package_dir/GeneSets.$long_name/data/ITFP.txt";
    open(IN, $input_file) || die "Could not open $input_file: $!\n";
    open(OUT, ">$output_file") || die "Could not write $output_file: $!\n";
    print OUT join("\t", qw(geneID termID termName dbName description))."\n";
    while(<IN>) {
        chomp;
        my ($tf, $target) = split(/\t/);
        exists($bad_symbols{$target}) && next;
        if (!exists($$symbol2entrezid{$target})) {
            #print STDERR "Unknown symbol $target\n";
            $bad_symbols{$target} = -1;
            next;
        } elsif (!defined($$symbol2entrezid{$target})) {
            #print STDERR "Ambiguous symbol $target\n";
            $bad_symbols{$target} = -1;
            next;
        }
        my $gene_id = $$symbol2entrezid{$target};
        my $term_id = "ITFP:$tf";
        my $name = "$tf targets";
        my $description = "Predicted targets of the ".
                $$tf_data{$tf}{"protein_symbol"}.
                " ($tf) transcription factor.";
        print OUT join("\t", $gene_id, $term_id, $name, "ITFP",
                $description)."\n";
    }
    close IN;
    close OUT;
}

#===============================================================================
# Main
#===============================================================================
my @organisms = &read_organism_file($ARGV[0]);
my $package_dir = $ARGV[1];
my %name_map = ("Homo sapiens" => "human", "Mus musculus" => "mouse",
        "Rattus norvegicus" => "rat");
foreach my $organism_name (@organisms) {
    my %symbol2entrezid = &create_symbol2id($organism_name);
    my %tf_data = &read_tf_data($name_map{$organism_name});
    &process_target_file($organism_name, $name_map{$organism_name},
            $package_dir, \%tf_data, \%symbol2entrezid);
}
