#!usr/bin/perl
use warnings ;
use strict ;

my $InputFile = $ARGV[0] ;
my $InputTable = $ARGV[1] ;
my $min_MAF = $ARGV[2] ;
my $NumWindows = $ARGV[3] ;
my $TimeOrientation = $ARGV[4] ;

#### MAKE SURE INPUT PARAMETERS SPECIFIED, INPUT FILES EXIST
unless(scalar @ARGV == 5){
	print "You didn't specify all four arguments the program needs:\n" ;
	print "\t 1.) Filename of alignment (in fasta format)\n" ;
	print "\t 2.) Filename of the tab-delimited table that specifies the sample date for each individual\n" ;
	print "\t 3.) Filename of alignment (in fasta format)\n" ;
	print "\t 4.) The number of windows, which should be less than half the sequence length or mayhem will ensue\n" ;
	print "\t 5.) Whether the program looks for the decay of linkage forward-in-time (F) or reverse-in-time (R),\n" ;
	exit ;
}
unless(-e "${InputFile}"){
	print "Couldn't detect input fasta file (1st argument), make sure the name is correct\n" ;
	exit ;
}
unless(-e "${InputTable}"){
	print "Couldn't detect input table with collection dates (2nd argument), make sure the name is correct\n" ;
	exit ;
}
unless( $min_MAF >= 0 && $min_MAF <= 0.5 ){
	print "The minor allele frequency cutoff (3rd argument) is not between 0 and 0.5\n" ;
	exit ;
}
unless($TimeOrientation eq "F" || $TimeOrientation eq "R" ){
	print "The time orientation parameter (5th argument) wasn't specified as \"F\" or \"R\", the only two options\n" ;
	exit ;
}

#############################################
# DATA STRUCTURES
#############################################
my $num_sites ;
my $sum ;
my %Seq_data ; #  $Seq_data{index} = scalar of 0's and 1's
my %Dates ;
my %Segsites_BiAllelic ; # $Segsites{index} = site ;
my %SumStat_Results ;
my %PairWise_PC ;

#############################################
# ASSEMBLE SEQUENCES PER IND
#############################################
open IN, "<${InputFile}" ;
my $CurrentInd = "" ;
while(<IN>){
	$_ =~ s/\r//g ;
	$_ =~ s/\n//g ;
	chomp $_ ;
	 
	if($_ =~ m/^>/){
		$_ =~ s/>//g ;
		$CurrentInd = $_ ;
		next ;
	}else{
		if( !exists $Seq_data{ $CurrentInd } ){
			$Seq_data{ $CurrentInd } = $_ ;
		}else{
			$Seq_data{ $CurrentInd } = $Seq_data{ $CurrentInd }.${_} ;
		}
	}
}
#############################################
# GET DATES FOR EACH INDIVIDUAL
#############################################
open IN, "<${InputTable}" ;
while(<IN>){
	$_ =~ s/\r//g ;
	$_ =~ s/\n//g ;
	chomp $_ ;
	
	my @line = split(/\t/, $_) ;
	$Dates{ $line[0] } = $line[1] ;
}

##### MAKE SURE INFORMATION FROM INPUT FILES CHECKS OUT
my $NumInd_InFile = scalar keys %Seq_data ; 
print "I counted ${NumInd_InFile} sequences in your input file, I hope this is correct!\n" ;
my $NumInd_InTable = scalar keys %Dates ; 
print "I counted ${NumInd_InTable} sequences with collection dates in your input table, this should be the same as the number of individuals\n" ;
unless($NumInd_InFile == $NumInd_InTable){
	print "The number of sequences in your fasta file and your table of collection dates do not match!\n" ;
	print "Please fix this and try again!" ;
	exit ;
}

foreach my $ind (sort {$a cmp $b} keys %Seq_data){
	if( !exists $Dates{$ind} ){
		print "I couldn't find a date for one or more sequences in your fasta file, \n" ;
		print "\tI don't know exactly why but make sure the names are the same in your fasta file and collection dates table\n" ;
		exit ;
	}
}
my %SeqLengths ;
foreach my $ind (keys %Seq_data){
	$num_sites = length $Seq_data{$ind} ;
	$SeqLengths{$num_sites}++ ;
}
if(scalar keys %SeqLengths >1){
	print "Not all of your sequences are the same length!\n" ;
	print "Fix this and try again.\n" ;
	exit ;
}

#############################################
## CONSTRUCT SINGLE SITE %AFS FOR SUMSTATS AND LD
#############################################
my $seg_site_index = 0 ;
foreach my $site ( 0 .. ($num_sites-1) ){
	my $af = 0 ;
	my %poly_bases ;
    my $indel = 0 ;
	foreach my $ind (keys %Seq_data){
		my $base = substr($Seq_data{$ind}, $site, 1) ;
		if( $base !~ m/-/ || $base !~ m/N/ ){
			$base = uc($base) ; #capitalize: codon table all caps
			$poly_bases{$base} ++ ;
		}else{
			$indel ++ ;
		}
	}
	# Is the site bi-allelic and not have any missing data?
	if( (scalar keys %poly_bases == 2) && ($indel==0) ){
		# Given it's polymorphic, does it pass the allele frequency cutoff specified?
		(my $first) = sort keys %poly_bases ;
		my $AF = $poly_bases{$first}/$NumInd_InFile ;
		if( $AF >= $min_MAF && $AF <= (1-$min_MAF) ){
			$Segsites_BiAllelic{$seg_site_index} = $site  ;
			$seg_site_index++ ;
		}		
	}
}

print "Number of sites: ", $num_sites, "\n" ;
print "Number of polymorphic sites: ", scalar keys %Segsites_BiAllelic, "\n" ;


#############################################
## CONSTRUCT WINDOWS
#############################################
my $window_size = $num_sites/$NumWindows ;
my %Windows ; # @{$Windows{$window}} = (all the counts)
foreach my $cnt ( sort{$a <=> $b} keys %Segsites_BiAllelic ){
    my $site = $Segsites_BiAllelic{$cnt} ;
    my $win = int($site/$window_size) ;
    push @{$Windows{$win}}, $cnt  ;
}

#foreach my $win (keys %Windows){
#	print $win, "\t", "@{$Windows{$win}}", "\n" ;
#}

my %DataSetsToIterateThrough ;
my %AnchorInds ;
my $first3 = 0 ;
if($TimeOrientation eq "F"){
	foreach my $ind (sort{$Dates{$a} <=> $Dates{$b}} keys %Dates){
		$AnchorInds{$ind}++ ;
		$first3++ ;
		if($first3==3){
			last ;
		}
	}
}elsif($TimeOrientation eq "R"){
	foreach my $ind (sort{$Dates{$b} <=> $Dates{$a}} keys %Dates){
		$AnchorInds{$ind}++ ;
		$first3++ ;
		if($first3==3){
			last ;
		}
	}

}
my $counter = 1 ;
foreach my $ind (sort {$a cmp $b} keys %Seq_data){
	if( !exists $AnchorInds{$ind} ){
		foreach my $x (sort {$a cmp $b} keys %AnchorInds){
			push @{$DataSetsToIterateThrough{$counter}}, $x ;
		}
		push @{$DataSetsToIterateThrough{$counter}}, $ind ;
		$counter++ ;
	}
}
print "These three sequences (below) were used to test whether a fourth sequence exhibited a 4th hapotype for each SNP pair:\n" ;
foreach my $ind (sort {$a cmp $b} keys %AnchorInds){
	print "\t", $ind, "\n" ;
}


#foreach my $key (sort{$a<=>$b} keys %DataSetsToIterateThrough){
#	print $key, "\t", "@{$DataSetsToIterateThrough{$key}}", "\n" ;
#}


#############################################
## CONSTRUCT %Pairwise_LD_LongRange, %Pairwise_Compatibility_LongRange values used for all subsequent calculations
#############################################


my ($lower, $upper) ;
foreach my $win (sort{$a <=> $b} keys %Windows){
	my @counts = @{$Windows{$win}} ; # these are $Biallelic_Syn_SNP_Ref_Pos count indices
	for ($lower = 0; $lower<((scalar @counts)-1); $lower++){
		my $site1 = $Segsites_BiAllelic{$lower} ;
		for ($upper = ($lower+1); $upper<(scalar @counts); $upper++){
			my $site2 = $Segsites_BiAllelic{$upper} ; 
			my %haplos ;
			foreach my $ind (keys %AnchorInds){
				my $b1 = uc(substr($Seq_data{$ind}, $site1, 1)) ;
				my $b2 = uc(substr($Seq_data{$ind}, $site2, 1)) ;
				$haplos{$b1.$b2}++ ; 		
			}
			# only go through all remaining individuals not originally included in the trio
			# if at least 3 haplotypes are seen in the "anchor" individuals, otherwise you'll never observe
			# a fourth haplotype!
			if(scalar keys %haplos < 3){
				# go through individuals of lower samp size gene, if they differ  
				foreach my $dataset ( sort {$a<=>$b} keys %DataSetsToIterateThrough ){
					my $Ind4 = ${$DataSetsToIterateThrough{$dataset}}[3] ;
					#print $Ind4, "\n" ;
					my $b1 = uc(substr($Seq_data{$Ind4}, $site1, 1)) ;
					my $b2 = uc(substr($Seq_data{$Ind4}, $site2, 1)) ;					
					$haplos{$b1.$b2}++ ; 		
					my $compatibility ;
					if((scalar keys %haplos) == 4){
						$compatibility = 0 ;
					}else{
						$compatibility = 1 ;
					}
		
					push @{$PairWise_PC{$win}{$dataset}}, $compatibility ;
				
					if( $haplos{$b1.$b2} > 1 ){
						$haplos{$b1.$b2}-- ;
					}else{
						delete $haplos{$b1.$b2} ;
					}
				}
			}else{
				next ;
			}
		}
	}
}


my %PC_mean_results ;

open OUT, ">./PCvsDate.txt" ;
print OUT "Window", "\t", "Sequence", "\t", "Date", "\t", "MeanPC", "\n" ;
foreach my $win (sort{$a<=>$b} keys %PairWise_PC){
	foreach my $dataset (sort{$a<=>$b} keys %{$PairWise_PC{$win}} ){
		print OUT $win, "\t", ${$DataSetsToIterateThrough{$dataset}}[3], "\t", $Dates{ ${$DataSetsToIterateThrough{$dataset}}[3] },  "\t" ;
		$sum = 0 ;
		foreach (@{$PairWise_PC{$win}{$dataset}}){
			$sum += $_ ;
		}
		print OUT $sum/scalar(@{$PairWise_PC{$win}{$dataset}}), "\n" ;	
	}
}
close OUT ;

exit ;

