use strict;
use Khanal;


my $start = time();

#check to see if correct number of command-line arguments have been given
unless (scalar @ARGV == 3)
{
	print "Usage: perl Khanal_Cliques_Bipartite.pl inputfile.txt output_cliques.txt output_bipartites.txt\n";
	exit;
}

#store command-line arguments in variables
my $input = $ARGV[0];
my $output_cliques = $ARGV[1];
my $output_bipartites = $ARGV[2];
chomp $output_bipartites;

#################################################################################
#create variables to store different data

#hash to store PPIs
my %interactions;

#hash to store proteins
my %hproteins;

#array to store proteins
my @proteins;

#variable to store the number of PPIs
my $ppicount = 0;

#array to store the cliques 
my @cliques;

#array to store the bipartites
my @bipartites;

##################################################################################

#load the data from input file
loadData($input, \%interactions, \%hproteins, \@proteins, \$ppicount);


##################################################################################

#Find the cliques of size-4 in a given network
findCliquesOfSize4(\%interactions, \@proteins, \@cliques, $output_cliques);

##################################################################################

#Find the cliques which form a output_bipartites
findBipartites(\@cliques, \%interactions,\@bipartites, $output_bipartites);

##################################################################################

my $end = time();
my $timeTaken = $end - $start;
print "The total time taken is: \t$timeTaken seconds\n";
