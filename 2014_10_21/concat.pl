#!/usr/bin/perl
use warnings;
#NOTES
#---------------------------------------------------------------------------------------#
## Developed by George Tiley - gtiley@ufl.edu
## usage: ./concat.pl .align.fasta (or however you end the files you want to align
#This script assumes alignments only differ by some prefix (e.g. 1.align.fasta, 2.align.fasta, 3.align.fasta ...) 
#This script reads through multiple aligned fasta files and concatenates them
#Missing data is inserted for missing taxa in each alignment
#All files need to be in the same folder as well as this script
#---------------------------------------------------------------------------------------#
#USER VARIABLES
#---------------------------------------------------------------------------------------#
#Define what character you would like to insert for missing data
$insert = "N";
#Name the output fasta file
$concat_file = "Concat.fasta";
#Name the log file, this tracks the number of taxa in each alignment and missing data
$concat_log = "Concat_log.txt";
#Get a phylip formatted file and a partition file for raxml
$phy_file = "Concat.phy";
$partition_file = "Concat.partitions.txt";
$type = "DNA"; #Can be DNA or AA
$Codon = 1; #Can be yes = 1  or no = 0
$AA_Q = "LGX"; #AA models need to be specified for each partition with RAxML. This assumes the same model for each partition.
#Write a nexus file with the MB block - you will want to finish the lset and mcmc commands
$MB_file = "Concat.MB.nex";
#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#
$SUFFIX =  "$ARGV[0]";
system "ls *.$SUFFIX > Alignments";
open FH1, "<Alignments";
open OUT1, ">$concat_file";
open OUT2, ">$concat_log";
open OUT3, ">$partition_file";
open OUT4, ">$phy_file";
open OUT5, ">$MB_file";
print OUT2 "TAXON\tMISSING_GENES\t\%_MISSING_GENES\t\%_MISSING_DATA\n";
%alltax = ();
%allgenes = ();
%allhash = ();
%misshash = ();
@gene_array = ();
$total_len = 0;
$N_tax = 0;
while (<FH1>)
{
	if (/((\S+)\.$SUFFIX)/)
	{
		$file = $1;
		$gene = $2;
		push @gene_array, $gene;
		if (! exists $allgenes{$gene})
		{
			$allgenes{$gene} = 1;
		}
		open FH2, "$file";
		while (<FH2>)
		{
			if (/\>(\S+)/)
			{
				$taxon = $1;
				if (! exists $alltax{$taxon})
				{
					$alltax{$taxon} = 1;
				}
			}
			elsif (/(\S+)/)
			{
				$allhash{$taxon}{$gene} = $1;
				if (! exists $misshash{$gene})
				{
					$misshash{$gene} = length ($allhash{$taxon}{$gene});
				}
				elsif (exists $misshash{$gene})
				{
					if ($misshash{$gene} < length ($allhash{$taxon}{$gene}))
					{
						$misshash{$gene} = length ($allhash{$taxon}{$gene});
					}
				}
			}
		}
		$total_len = $total_len + $misshash{$gene};
		close FH2,
	}
}
close FH1;
$Ngenes = scalar @gene_array;
%missing_data = ();
foreach $TAX (keys %alltax)
{
	print OUT1 ">$TAX\n";
	print OUT2 "$TAX\t";
	$miss_gene_count = 0;
	$miss_data_count = 0;
	$N_tax++;
	foreach $GENE (keys %allgenes)
	{
		if (! exists $allhash{$TAX}{$GENE})
		{
			$miss_gene_count++;
			for (0..(($misshash{$GENE})-1))
			{
				print OUT1 "$insert";
				$miss_data_count++;
			}
		}
		if (exists $allhash{$TAX}{$GENE})
		{
			print OUT1 "$allhash{$TAX}{$GENE}";
		}
	}
	print OUT1 "\n";
	$missing_genes = ($miss_gene_count / $Ngenes) * 100;
	$missing_data = ($miss_data_count / $total_len) * 100;
	print OUT2 "$miss_gene_count\t$missing_genes\t$missing_data\n";
}	
close OUT1;
close OUT2;
#----------------------------------------------------------------------------------------#
open FH2, "<$concat_file";
print OUT4 "$N_tax\t$total_len\n";
print OUT5 "#NEXUS\nbegin data;\ndimensions ntax=$N_tax nchar=$total_len;\nformat missing=$insert gap=- matchchar=. datatype=$type;\nmatrix\n";
while (<FH2>)
{
	if (/>(\S+)/)
	{
		$TAXON = $1;
	}
	elsif (/(\S+)/)
	{
		$SEQ = $1;
		print OUT4 "$TAXON\t\t$SEQ\n";
		print OUT5 "$TAXON\t\t$SEQ\n";
	}
}	
close FH2;
print OUT5 ";\nend;\n";
#----------------------------------------------------------------------------------------#
$start = 0;
print OUT5 "BEGIN mrbayes;\n";
@Partitions = ();
if ($type eq "DNA")
{
	if ($Codon == 0)
	{
		foreach $gene (keys %misshash)
		{
			$begin = ($start + 1);
			$end = ($begin + ($misshash{$gene} - 1));
			print OUT3 "$type, $gene= $begin-$end\n";
			print OUT5 "charset $gene = $begin-$end;\n";
			push @Partitions, $gene;
			if ($start < ($total_len + 1))
			{
				$start = $end;
			}
		}
		$N_Partitions = scalar @Partitions;
		$last = $N_Partitions - 1;
		print OUT5 "partition Names = $N_Partitions: ";
		for (0..($N_Partitions - 2))
		{
			$i = $_;
			print OUT5 "$Partitions[$i], ";
		}
		print OUT5 "$Partitions[$last];\n";
		print OUT5 "set partition = Names;\n"; 
		print OUT5 "lset applyto=(";
		for (1..($N_Partitions - 1))
		{
			$i = $_;
			print OUT5 "$i,";
		}
		print OUT5 "$N_Partitions) nst=6 rates=gamma;\n";
		print OUT5 "mcmc ngen=10000000 nruns=4 nchains=4 samplefreq=100;\nend;";
	}
	if ($Codon == 1)
	{
		foreach $gene (keys %misshash)
		{
			$begin = ($start + 1);
			$pos2 = $begin + 1;
			$pos3 = $begin + 2;
			$end = ($begin + ($misshash{$gene} - 1));
			$POS1 = $gene . "_Pos1";
			$POS2 = $gene . "_Pos2";
			$POS3 = $gene . "_Pos3";
			print OUT3 "$type, $POS1= $begin-$end\\3\n";
			print OUT3 "$type, $POS2= $pos2-$end\\3\n";
			print OUT3 "$type, $POS3= $pos3-$end\\3\n";
			print OUT5 "charset $POS1 = $begin-$end\\3;\n";
			print OUT5 "charset $POS2 = $pos2-$end\\3;\n";
			print OUT5 "charset $POS3 = $pos3-$end\\3;\n";
			push @Partitions, "$POS1";
			push @Partitions, "$POS2";
			push @Partitions, "$POS3";
			if ($start < ($total_len + 1))
			{
				$start = $end;
			}
		}
		$N_Partitions = scalar @Partitions;
		$last = $N_Partitions - 1;
		print OUT5 "partition Names = $N_Partitions: ";
		for (0..($N_Partitions - 2))
		{
			$i = $_;
			print OUT5 "$Partitions[$i], ";
		}
		print OUT5 "$Partitions[$last];\n";
		print OUT5 "set partition = Names;\n"; 
		print OUT5 "lset applyto=(";
		for (1..($N_Partitions - 1))
		{
			$i = $_;
			print OUT5 "$i,";
		}
		print OUT5 "$N_Partitions) nst=6 rates=gamma;\n";
		print OUT5 "mcmc ngen=10000000 nruns=4 nchains=4 samplefreq=100;\nend;";
	}
}
if ($type eq "AA")
{
	foreach $gene (keys %misshash)
	{
		$begin = ($start + 1);
		$end = ($begin + ($misshash{$gene} - 1));
		print OUT3 "$AA_Q, $gene= $begin-$end\n";
		print OUT5 "charset $gene = $begin-$end;\n";
		push @Partitions, $gene;
		if ($start < ($total_len + 1))
		{
			$start = $end;
		}
	}
	$N_Partitions = scalar @Partitions;
	$last = $N_Partitions - 1;
	print OUT5 "partition Names = $N_Partitions: ";
	for (0..($N_Partitions - 2))
	{
		$i = $_;
		print OUT5 "$Partitions[$i], ";
	}
	print OUT5 "$Partitions[$last];\n";
	print OUT5 "set partition = Names;\n"; 
	print OUT5 "prset applyto=(";
	for (1..($N_Partitions - 1))
	{
		$i = $_;
		print OUT5 "$i,";
	}
	print OUT5 "$N_Partitions) aamodelpr=fixed(WAG);\n";
	print OUT5 "mcmc ngen=10000000 nruns=4 nchains=4 samplefreq=100;\nend;";
}
close OUT3;
exit;