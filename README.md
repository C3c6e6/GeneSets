# What does this package do?

This package is a set of scripts to build gene set annotation packages for different organism for use with the SetRank BioConductor package.

# How do I install it?

There is no installation. Just clone this repository to your local machine. Make sure all scripts in the bin folder have their execute permission set. These scripts need the following requirements:
* R 3.0.1 or higher with the [`GO.db`](http://bioconductor.org/packages/release/data/annotation/html/GO.db.html) Bioconductor package installed.
* Perl 5.12 or higher with the [`Date::Format`](http://search.cpan.org/~gbarr/TimeDate-2.30/lib/Date/Format.pm) module installed.
* The `wget` utility
* GNU Make

# How do I use it?

## Building the GeneSets packages

The only thing you have to do is tell the package for which organisms you which to build annotation tables. In root of the repository, there is a file called `organisms`. Open this file with your favorite text editor and add the names of the organism(s) you want to build annotation tables for and remove any names you're not interested in. Make sure there is only one name per line and that the lines do not contain any leading or trailing white space characters. Also, make sure that the names you are using correspond exactly to the official names used by the NCBI taxonomy database, as this is the reference used by the GeneSets package. Using the correct name is especially important when working with bacterial strains and substrains. For instance, if you want to create annotation tables for Escherichia coli K12, substrain MG1655, you should write exactly: `Escherichia coli str. K-12 substr. MG1655`. When in doubt, you can look up the name of your species on the NCBI website at http://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi and simply copy-and-paste the name into your editor.

Once you've edited the `organisms file`, you are ready to build the packages. All you have to do is open your UNIX-shell, go to the directory of the repository:
```
make
```

Issuing this command will start the build process, which might take a while to complete, mainly depending on your network connection speed. 

Once the process is finished you should see, for each organism listed in the organisms file, a file appear named `GeneSets.Genus.species_YY.MM.DD.tar.gz`, with `Genus.species` the name of the organism and `YY.MM.DD` the date the file, e.g. `GeneSets.Homo.sapiens_13.10.18.tar.gz`. Each of these files is an R package -- the date string functions as version number -- that is ready to install.

To install several packages at once, you can use the following bash command:

```bash
for f in *_YY.MM.DD.tar.gz; do R CMD INSTALL $f; done
```

Of course, you have to replace `YY.MM.DD` with the actual date string.

## Updating the annotation packages

To keep your collection of annotation table packages up-to-date, all you have to do is to regularly run the make command. The GeneSets package will automatically detect if the gene set annotation data for a given organism has changed compared to the last time you ran make and create an updated version if necessary. Don't forget to install the updated packages afterwards.


