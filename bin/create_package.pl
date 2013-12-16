#!/usr/bin/env perl

use strict;
use warnings;
use File::Spec;
use Date::Format;
use File::Copy;

my $template_dir = "templates/";
my $package_dir = $ARGV[0];
$package_dir =~ s/\/$//g;

#===============================================================================
# Subroutines
#===============================================================================
sub read_file {
    my $file_name = $_[0];
    open(IN, $file_name) || die "Could not open $file_name: $!\n";
    my @contents = <IN>;
    close IN;
    return join("", @contents);
}

sub get_object_list {
    my $data_dir = $_[0];
    my @object_list = <$data_dir/*.txt>;   
    foreach my $object (@object_list) {
        $object =~ s/.+\/(\w+)\.txt$/$1/;
    }
    return(@object_list);
}

sub write_documentation {
    my @objects = @_;
    my $documentation = &read_file("$template_dir/documentation.Rd");
    my $aliases = join("\n", map { "\\alias{$_}" } "allDBs", @objects);
    my $organism_name = (File::Spec->splitdir($package_dir))[-1];
    $organism_name =~ s/GeneSets_//;
    $organism_name =~ s/_/ /g;
    $documentation =~ s/__ORGANISM__/$organism_name/;
    $documentation =~ s/__ALIASES__/$aliases/;
    open(OUT, ">$package_dir/man/documentation.Rd") || 
            die "Could not write $package_dir/man/documentation.Rd: $!\n";
    print OUT $documentation;
    close OUT;
}

sub write_description {
    my ($description, $version,$date, $description_file, $package_name);
    $description = &read_file("$template_dir/DESCRIPTION");
    $description_file = "$package_dir/DESCRIPTION";
    $version = time2str("%y.%m.%d", time);
    $version =~ s/\b0//g;
    $date = time2str("%Y-%m-%d", time);
    $package_name = (File::Spec->splitdir($package_dir))[-1];
    $description =~ s/__PKGNAME__/$package_name/;
    $description =~ s/__VERSION__/$version/;
    $description =~ s/__DATE__/$date/;
    open(OUT, ">$description_file") ||
            die "Could not write $description_file: $!\n";
    print OUT $description;
    close OUT;
}

sub write_namespace {
    copy("$template_dir/NAMESPACE", "$package_dir/NAMESPACE")
            || die "Could not copy NAMESPACE file: $!\n";
}

sub create_R_objects {
    my @object_names = @_;
    my $out_dir = "$package_dir/data";
    my $arguments = 'header=TRUE, sep="\\t", quote="", colClasses="factor"';
    open(R, "|R --no-save --slave") || die "Could open a pipe to R: $!\n";
    foreach (@object_names) {
        print R "$_ = read.table(\"$out_dir/$_.txt\", $arguments)\n";
    }
    my $object_string = join(", ", @object_names);
    print R "allDBs = rbind($object_string)\n";
    print R "save(allDBs, $object_string, file=\"$out_dir/data.rda\")\n";
    close R;
    unlink( map {"$out_dir/$_.txt"} @object_names );
}

sub build_package {
    my $cur_dir = File::Spec->curdir();
    my @build_dir = File::Spec->splitdir($package_dir);
    my $package_name = pop(@build_dir);
    chdir(File::Spec->catdir(@build_dir));
    system("R CMD build $package_name");
    chdir($cur_dir);
}

#===============================================================================
# Main
#===============================================================================
my @object_list = &get_object_list("$package_dir/data/");

&write_documentation(@object_list);
&write_description();
&write_namespace();
&create_R_objects(@object_list);
&build_package();
