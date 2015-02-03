use strict;

##############################################################################
# Subroutine to load the protiens and the PPIs into memory
# first parameter: filename of the input file 
# second parameter: reference to a hash where PPIs are stored
# third parameter: reference to an hash where proteins are stored
# fourth parameter: reference to an array where proteins are stored
# fifth parameter: reference to a scalar to store number of PPIs
##############################################################################
sub loadData
{
	my $filename = $_[0];
	my $interaction = $_[1];
	my $proteins = $_[2];
	my $prot = $_[3];
	my $ppicount = $_[4];

	#open the input file
	unless (open(FHR, $filename))
	{
	print "Couldn't open input file";
	exit;
	}

	#input the data into different variables in different forms
	while (my $line = <FHR>)
	{
		my @data = split ("\t", $line);
		chomp @data;
		
		my @individual;
		push(@individual, $data[0]);
		push(@individual, $data[1]);

		@individual = sort @individual;

		#store the individual proteins in a hash
		$$proteins{$individual[0]}++;
		$$proteins{$individual[1]}++;

		#store the interaction or edges in a hash
		my $ppi = $individual[0]."\t".$individual[1];
		$$ppicount++;
		$$interaction{$ppi}++;
	}
	#store the proteins in an sorted array
	@$prot = keys %$proteins;
	@$prot = sort @$prot;

	#print some useful information
	print "The number of PPI is: \t $$ppicount \n";
	print "The number of proteins is:\t", scalar @$prot, "\n";
}


###############################################################################
#Subroutine to find the cliques of size 4 in a given PPIs
#first parameter: reference to hash of PPIs
#second parameter: reference to array of proteins
#third parameter: reference to array to store cliques
#fourth parameter: output filename 
###############################################################################

sub findCliquesOfSize4
{
	my $ppi = $_[0];
	my $prot = $_[1];
	my $cliques = $_[2];
	my $output = $_[3];

	#declare an array to store the edges
	my @edges;
	
	#for-loop to get the first node
	for(my $i = 0; $i < scalar @$prot - 3; $i++)
	{
		#push the first node into an array
		my @nodes;
		push (@nodes, $$prot[$i]);

		#for-loop to get the second node
		for(my $j = $i +1; $j < scalar @$prot - 2; $j++)
		{
			#push the second node
			push(@nodes, $$prot[$j]);

			#for-loop to get the third node
			for(my $k = $j + 1; $k < scalar @$prot - 1; $k++)
			{
				#push the third node
				push(@nodes, $$prot[$k]);

				#for-loop to get the fourth node
				for(my $l = $k + 1; $l < scalar @$prot; $l++)
				{
					#push the fourth node
					push(@nodes, $$prot[$l]);

					#check if @nodes form a clique of size 4
					my $check = isClique(\@nodes, $ppi, \@edges);
					if ($check == 1)
					{
						#push the combination of nodes into cliques array
						push(@$cliques, join("\t", @nodes));
					}
					pop @nodes;
				}
				pop @nodes; 
			}
			pop @nodes;
		}
		pop @nodes;
	}

	#print some useful information:
	print "The number of clique of size 4 is:\t", scalar @$cliques, "\n";

	#write output to a output file
	writeCliques($cliques, \@edges, $output);

}



######################################################################################
# Subroutine to check if a given array of nodes form a clique of size n given the PPIs
# first argument : reference to array of nodes to be checked
# second argument : referecne to hash of PPIs
# third argument : reference to an array to store all the edges
######################################################################################
sub isClique
{
	#store the given arguments
	my $nodes = $_[0];
	my $PPI = $_[1];
	my $edges = $_[2];

	my $allEdges = ""; 

	#for-loop to create all the possible edges
	for(my $k = 0; $k < scalar @$nodes - 1; $k++)
	{
		for(my $l = $k + 1; $l < scalar @$nodes; $l++)
		{
			
			my $edge = $$nodes[$k]."\t".$$nodes[$l];

			#check if the node exists in %PPI and return 0 if it does not
			unless (exists $$PPI{$edge})
			{
				return 0;
			}

			if($allEdges eq "")
			{ 
				$allEdges = $edge;
			}
			else
			{
				$allEdges = $allEdges.",".$edge;
			}

		}
	}
	push(@$edges, $allEdges);

	return 1;
}



#######################################################################################
# Subroutine to find all the bipartites
# first argument: reference to an array of cliques
# second argument: reference to a hash of PPI
# third argument: reference to an array to store bipartites
# fourth argument: outfile name to write bipartites and the edges
#######################################################################################
sub findBipartites
{
	my $cliques = $_[0];
	my $ppi = $_[1];
	my $bipartites = $_[2];
	my $outfile = $_[3];
	my $numBipartite = 0;

	#declare array to store all the edges to form a bipartite
	my @edges;
	#go through all the cliques to see if they form bipartite
	for(my $i = 0; $i < scalar @$cliques - 1; $i++)
	{
		my @clique1 = split ("\t", $$cliques[$i]);
		for(my $j = $i + 1; $j <scalar @$cliques; $j++)
		{
			
			my @clique2 = split ("\t", $$cliques[$j]);



			#check if they form a bipartite
			my $check = isBipartite(\@clique1, \@clique2, $ppi, \@edges);
			if ($check == 1)
			{
				my $bipartite = @$cliques[$i].":".@$cliques[$j];
				push(@$bipartites, $bipartite);
				$numBipartite++;
			}
		}
	}
	#print some important information
	print "The number of bipartites is: \t", $numBipartite, "\n";

	#print bipartite and edges in an output file
	writeBipartites($bipartites, \@edges, $outfile);
}


########################################################################
# Subroutine to check whether given two cliques form a bipartite
# first argument: ref. to clique1
# second argument: ref. to clique2
# third argument: ref. to hash of interactions
# fourth argument: ref. to store all the edges to form the bipartite
########################################################################
sub isBipartite
{
	my $clique1 = $_[0];
	my $clique2 = $_[1];
	my $PPI = $_[2];
	my $edges = $_[3];
	#string to store all the edges present to form the bipartite
	my $allEdges = "";

	#hash to store the count of edges for each protein
	my %nodes;

	#check if the nodes are the same
	for (my $a = 0; $a < scalar @$clique1; $a++)
	{
		$nodes{$$clique1[$a]} = 0;
		$nodes{$$clique2[$a]} = 0;
	}

	#for-loop to check all the combinations of cliques
	for(my $i = 0; $i < scalar @$clique1; $i++)
	{
		for(my $j = 0; $j < scalar @$clique2; $j++)
		{
			if($$clique1[$i] eq $$clique2[$j])
			{
				return 0;
			}

			else
			{
				
				my $edge = join("\t", sort($$clique1[$i], $$clique2[$j]));

				if (exists $$PPI{$edge})
				{

					$nodes{$$clique1[$i]}++;
					$nodes{$$clique2[$j]}++;

					if ($allEdges eq "")
					{
						$allEdges = $edge;
					}
					else
					{
						$allEdges = $allEdges.",".$edge;
					}
				}
				
			}
		}
	}


	#check if at least one connection exists
	my @numEdges = values %nodes;
	
	for (my $i = 0; $i < scalar @numEdges; $i++)
	{
		
		if ($numEdges[$i] == 0)
		{
			return 0;
		}
	}
	push(@$edges, $allEdges);
	return 1;
}



##################################################################################
# Subroutine to write the cliques and the edges to an output file
# first argument: reference to an array of cliques
# second argument: reference to an array of edges for cliques
# third argument: name of output file
##################################################################################
sub writeCliques
{
	#store the arguments
	my $cliques = $_[0];
	my $edges = $_[1];
	my $output = $_[2];

	#open the output file to write the cliques 
	unless (open(CFHR, ">$output"))
	{
		print "Coudn't open output file for cliques\n";
		exit;
	}

	#writing the cliques and edges to the output file
	for(my $i = 0; $i < scalar @$cliques; $i++)
	{
		print CFHR "$$cliques[$i]:";
		print CFHR "$$edges[$i]\n";
	}

}

###################################################################################
# Subroutine to write the bipartites and the edges to an output file
# first argument: ref. to an array of bipartites
# second argument: ref. to an array of edges to form bipartite
# third argument: output file name
###################################################################################
sub writeBipartites
{
	#store the arguments
	my $bipartites = $_[0];
	my $edges = $_[1];
	my $output = $_[2];

	#open the output file to write the cliques 
	unless (open(BFHR, ">$output"))
	{
		print "Coudn't open output file for bipartites\n";
		exit;
	}

	#write the bipartites and edges to the output file
	for(my $i = 0; $i < scalar @$bipartites; $i++)
	{
		print BFHR "$$bipartites[$i]:";
		print BFHR "$$edges[$i]\n";
	}

}

1;
