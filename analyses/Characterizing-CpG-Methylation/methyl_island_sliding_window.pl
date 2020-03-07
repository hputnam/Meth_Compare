#!/usr/bin/perl -w

# usage: ./meth_island_sliding_window.pl [window size] [mCpG fraction] [sorted mCpG File]
# input: file with positions of mCpG in sorted order (col1 = scaffold, col2 = position of bp)
# output: prints the location of predicted methylated islands and number of mCpG's per island to STDOUT in following format: (col1 = scaffold, col2 = start position of island, col3 = end position of island, col4 = Number of mCpG's in the island)

use strict;
use warnings;

# Initialization
my $winSize = $ARGV[0];
my $mcpgFrac = $ARGV[1];
my $stepSize = $ARGV[2];

# Runs main subroutine
Main();

# Main sub 
sub Main {
	my @meth_file = read_in();
	my %meth_hash = hash_in(\@meth_file);
	sliding_win(\%meth_hash);
}

# read input file from command line into an array
sub read_in {
	die "$ARGV[3] cannot be found\n" unless (open(DATA, $ARGV[3])); 
	my @input1 = <DATA>;
	close (DATA);
	return @input1;
}

# Create hash of input file with key being scaffold name and a sorted array for each scaffold indexing mCpG position
sub hash_in {
	my @meth_file = @{$_[0]};
	
	my %meth_hash;
	my $k = 0;
	my $scaf = "";
	
	foreach my $i (@meth_file) {
		if ($i =~ /^(.*?)\t(.*)\n$/) {
			if ($1 eq $scaf) {$k++;}
				else {$k = 0;}
			$meth_hash{$1}[$k] = $2;
			$scaf = $1;
		}
	}
	return %meth_hash;
}

# Runs sliding window starting from an mCpG site. if it passes mCpG fraction threshold, it extends by the window size until threshold is no longer passed.
sub sliding_win {
	my %meth_hash = %{$_[0]};

	my @keys = sort keys %meth_hash;

	# Goes through each scaffold at a time
	foreach my $i (@keys) {

		for (my $j = 0; $j < scalar(@{$meth_hash{$i}}); $j++) {

			my $winStart = $meth_hash{$i}[$j]; # Start of region
			my $initialWinStart = $meth_hash{$i}[$j];
			my $winStop = $winStart + $winSize - 1; # End of region
			my $methNum = 1; # Number of methylated sites
			my $region = $winSize; # This variable will be extended if region meets criteria
			my $oldMethNum; # Stores number of mCG's from previous iteration of region
			my $oldj = $j; # Stores starting index

			EXTEND:	

			# Counts number of methylated CpG's in region	
			while (exists $meth_hash{$i}[$j] && $meth_hash{$i}[$j] >= $winStart && $meth_hash{$i}[$j] <= $winStop) {
				$methNum++;
				$j++;
			}
			$j--;
			$methNum--;

			# If region meets criteria by passing mCG fraction threshold, the region is extended and sent back to EXTEND to recount mCG's
			if ($methNum / $region >= $mcpgFrac && exists $meth_hash{$i}[$j + 1]) {
				$winStart = $initialWinStart + $region;
				$winStop = $winStart + $stepSize - 1;
				$region += $stepSize;
				$oldMethNum = $methNum; # Stores number of mCG's in the current iteration of the region before extending
				$oldj = $j; # Stores index of last mCG of the current iteration of the region before extending
				$j++;
				$methNum++;
				goto EXTEND;
				
				# Region is rejected			
			}	elsif ($methNum / $region < $mcpgFrac && $region == $winSize) {$j = $oldj; next;}

				# Region is accepted and has reached end of scaffold; Prints results
				elsif ($methNum / $region >= $mcpgFrac && ! exists $meth_hash{$i}[$j + 1]) {
					$region -= $stepSize;
					$winStart = $initialWinStart;
					$winStop = $meth_hash{$i}[$j]; # End site is stored as the last mCG in the region
					print "$i\t$winStart\t$winStop\t$methNum\n";
				}
				
				# Current region is rejected so reverts back to previous region iteration and prints results
				else {
					$region -= $stepSize;
					$winStart = $initialWinStart; 
					$winStop = $meth_hash{$i}[$oldj]; # End site is stored as the last mCG in the region
					$j = $oldj;
					print "$i\t$winStart\t$winStop\t$oldMethNum\n";
				}
		}
	}
}
