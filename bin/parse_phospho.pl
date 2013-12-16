#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use FileHandle;
my $db_name = "PhosphoSitePlus";

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

sub map_uniprot_ids {
    my @organism_names = @_;
    my $script = dirname($0)."/create_conversion_table.R";
    my (%uniprot2organism, %uniprot2entrezid);
    for my $organism (@organism_names) {
        my $annot_package = $organism;
        $annot_package =~ s/\s/./g;
        open(IN, "$script $annot_package UNIPROT ENTREZID 2> /dev/null |")
                || die "Could not open pipe: $!\n";
        <IN>; #Skip first line with column names.
        while (<IN>) {
            chomp;
            my ($uniprot, $entrezid) = split(/\t/);
            $uniprot2entrezid{$uniprot}=$entrezid;
            $uniprot2organism{$uniprot}=$organism;
        }        
        close IN;
    }
    return(\%uniprot2organism, \%uniprot2entrezid);
}

sub get_kinase_names {
    my @organism_names = @_;
    my $script = dirname($0)."/create_conversion_table.R";
    my %uniprot2name;
    for my $organism (@organism_names) {
        my $annot_package = $organism;
        $annot_package =~ s/\s/./g;
        open(IN, "$script $annot_package UNIPROT GENENAME 2> /dev/null |")
                || die "Could not open pipe: $!\n";
        <IN>; #Skip first line with column names.
        while (<IN>) {
            chomp;
            my ($uniprot, $gene_name) = split(/\t/);
            $uniprot2name{$uniprot}=$gene_name;
        }        
        close IN;
    }
    return \%uniprot2name;
}

sub process_kinase_data {
    my ($uniprot2organism, $uniprot2entrezid, $uniprot2name, $data_file,
           $read_flag, $package_dir);
    my %file_handles;
    ($uniprot2organism, $uniprot2entrezid, $uniprot2name, $package_dir) = @_;
    $data_file = "Kinase_Substrate_Dataset";
    $read_flag = 0;
    open(IN, $data_file) || die "Could not open $data_file: $!\n";
    while(<IN>) {
        if (m/^KINASE/) {
            $read_flag=-1;
            next;
        }
        $read_flag || next;
        chomp;
        my ($kinase_symbol, $kinase_acc, undef, undef, undef, undef, undef,
                $substrate_acc) = split(/\t/);
        exists($$uniprot2organism{$kinase_acc}) || next;
        exists($$uniprot2organism{$substrate_acc}) || next;
        my $organism = $$uniprot2organism{$kinase_acc};
        ($organism eq $$uniprot2organism{$substrate_acc}) || next;
        exists($file_handles{$organism}) || ($file_handles{$organism} =
                &create_file_handle($organism, $package_dir));
        my $gene_id = $$uniprot2entrezid{$substrate_acc};
        my $term_id = "$db_name:$kinase_symbol";
        my $term_name = "$kinase_symbol substrates";
        my $kinase_name = $$uniprot2name{$kinase_acc};
        my $description =
                "All substrates of the $kinase_name ($kinase_symbol).";
        $file_handles{$organism}->print(join("\t", $gene_id, $term_id,
                $term_name, $db_name, $description)."\n");
    }
    close IN;
    for my $fh (values %file_handles) {
        $fh->close();
    }
}

sub create_file_handle {
    my ($organism, $output_file, $output_dir, $fh);
    ($organism, $output_dir) = @_;
    $organism =~ s/\s/./g;
    $output_file = "$output_dir/GeneSets.$organism/data/$db_name.txt";
    $fh = FileHandle->new(">$output_file") ||
            die "Could not open $output_file: $!\n";
    my @column_names = qw(geneID termID termName dbName description);
    print $fh join("\t", @column_names)."\n";
    return($fh);
}

#===============================================================================
# Main
#===============================================================================
my ($uniprot2organism, $uniprot2entrezid, $uniprot2name);
my @organisms = &read_organism_file($ARGV[0]);
my $package_dir = $ARGV[1];
($uniprot2organism, $uniprot2entrezid) = &map_uniprot_ids(@organisms);
$uniprot2name = &get_kinase_names(@organisms);
&process_kinase_data($uniprot2organism, $uniprot2entrezid, $uniprot2name, 
        $package_dir);
