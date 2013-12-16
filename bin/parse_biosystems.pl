#!/usr/bin/env perl

use warnings;
use strict;
use FileHandle;

my $biosys_info_file = "bsid2info.gz";
my $biosys_gene_file = "biosystems_gene.gz";
my $taxonomy_file = "names.dmp";
my $gene_info_file = "gene_info.gz";
my $output_dir = "../../package_data/";

#===============================================================================
# Subroutines
#===============================================================================
sub read_organism_list {
    my $input_file = $_[0];
    open(IN, $input_file) || die "Could not open $input_file: $!\n";
    my @organisms = <IN>;
    close IN;
    chomp(@organisms);
    return(@organisms)
}

sub get_taxids {
    my @organism_list = @_;
    my (%organisms_whitelist, %organisms);
    @organisms_whitelist{@organism_list} = (-1) x scalar(@organism_list);
    open(IN, $taxonomy_file) || die "Could not open $taxonomy_file: $!\n";
    while (<IN>) {
        chomp;
        my ($id, $name, $other, $type) = split(/\t\|\t/);
        $type =~ s/\t\|$//;
        ($type eq "scientific name") || next;
        exists($organisms_whitelist{$name}) && ($organisms{$id} = $name);
    }
    close IN;
    return %organisms;
}

sub read_gene_ids {
    my %organisms = %{$_[0]};
    my %gene2tax_id;
    open(IN, "gunzip -c $gene_info_file|")
            || die "Could not open $gene_info_file: $!\n";
    while(<IN>) {
        chomp;
        my ($tax_id, $gene_id) = split(/\t/);
        exists($organisms{$tax_id}) || next;
        $gene2tax_id{$gene_id} = $tax_id;
    }
    close IN;
    return %gene2tax_id;
}

sub read_set_definitions {
    my %gene_sets = %{$_[0]};
    my %set_definitions;
    open(IN, "gunzip -c $biosys_info_file|") || 
            die "Could not open $biosys_info_file: $!\n";
    while(<IN>) {
        chomp;
        my %record;
        @record{qw(bsid db db_id name type scope tax_id description)} = 
                split(/\t/);
        exists($gene_sets{$record{"bsid"}}) || next;
        defined($record{"description"}) || ($record{"description"} = "");
        $set_definitions{$record{"db"}}{$record{"bsid"}} = \%record;
    }
    close IN;
    return %set_definitions;
}

sub read_set_memberships {
    my (%gene2tax_id, %gene_sets);
    %gene2tax_id = %{$_[0]};
    open(IN, "gunzip -c $biosys_gene_file|") || 
            die "Could not open $biosys_gene_file: $!\n";
    while (<IN>) {
        chomp;
        my ($bsid, $gene_id) = split(/\t/);
        exists($gene2tax_id{$gene_id}) || next;
        push @{$gene_sets{$bsid}}, $gene_id;
    }
    close IN;
    return %gene_sets;
}

sub create_file_handle {
    my ($db_name, $organism) = @_;
    $db_name =~ s/\s/_/g;
    $organism =~ s/\s/./g;
    my $output_file = "$output_dir/GeneSets.$organism/data/$db_name.txt";
    my $fh = FileHandle->new(">$output_file") ||
            die "Could not open $output_file: $!\n";
    my @column_names = qw(geneID termID termName dbName description);
    print $fh join("\t", @column_names)."\n";
    return($fh);
}


#===============================================================================
# Main
#===============================================================================
print STDERR "Reading taxonomy...";
my %organisms = &get_taxids(&read_organism_list($ARGV[0]));
print STDERR "\nReading gene IDs...";
my %gene2tax_id = &read_gene_ids(\%organisms);
print STDERR "\nReading set memberships...";
my %gene_sets = &read_set_memberships(\%gene2tax_id);
print STDERR "\nReading set definitions...";
my %set_definitions = &read_set_definitions(\%gene_sets);

print STDERR "\nWriting data...";
foreach my $db (keys %set_definitions) {
    my %file_handles;
    foreach my $bsid (keys %{$set_definitions{$db}}) {
        my %record = %{$set_definitions{$db}{$bsid}};
        my $term_id = $record{'db_id'};
        foreach my $gene_id (@{$gene_sets{$bsid}}) {
            my $organism = $organisms{$gene2tax_id{$gene_id}};
            exists($file_handles{$organism}) ||($file_handles{$organism} =
                    &create_file_handle($db, $organism)); 
            my @out = ($gene_id, $term_id, $record{"name"}, $record{"db"},
                  $record{"description"});
            $file_handles{$organism}->print(join("\t", @out)."\n");
        }
    }
    foreach my $fh (values %file_handles) {$fh->close()}
}
print STDERR "\n";
