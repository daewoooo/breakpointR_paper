#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use List::Util qw/shuffle/;
use constant PI=>3.14159265;

#my $GENOME = '/home/daewoooo/GENOMES/human_GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna';

##usage##
if(@ARGV < 1) {
	print "WRONG argument submitted: try simulateStrandSreads.pl -h\|-help\n";
	exit;
} elsif ($ARGV[0] eq '-h' || $ARGV[0] eq '-help') {
	help();
	exit;
}

my $fasta = $ARGV[0];
my $fraction = $ARGV[1];
my $background = $ARGV[2];
my $ID = $ARGV[3];

open my ($read1), '>', 'read1.fq' or die "Can't open file: $!\n";
open OUT, '>', $ID.'_SCE_list.txt' or die "Can't open file: $!\n";
print OUT join "\t", qw(Chr Pos Filename),"\n";

my $db = Bio::DB::Fasta->new($fasta);
my @ids = $db->get_all_ids;

my @chrom = ('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX');

my $genome_len = 0;

my $mean_frg_length = 150;
my $SD_frg_length = 10;
my $read_len = 76;

## Postions of centromes in acrocentromeric chromosomes
my %acrocent_chroms = (	'chr13' => [16000000, 19000000],
			'chr14' => [16000000, 19000000],
			'chr15' => [17000000, 20000000],
			'chr21' => [11288129, 14288129],
			'chr22' => [13000000, 16000000]
		      );

my $NumOfLibs = 10;
my %seen;
my $counter = 0;
my $genome = 0;

my $created_reads = 0;

foreach my $num (1..$NumOfLibs) {

	open my ($read1), '>', 'read1.fq' or die "Can't open file: $!\n";

	my $prefix = $ID;
	$prefix = $prefix.'_'.$num;
	warn "Working on $prefix ...\n";

	my @chrom_copy = @chrom;
	@chrom_copy = shuffle(@chrom_copy);
	my @wc = splice(@chrom_copy, 12);
	@chrom_copy = shuffle(@chrom_copy);
	my @w = splice(@chrom_copy, 6);
	my @c = @chrom_copy;
		
	@chrom_copy = @chrom;
	@chrom_copy = grep { $_ !~ /Y/ } @chrom_copy;
	@chrom_copy = shuffle(@chrom_copy);

	my @SCE = splice(@chrom_copy, 0, 5);

foreach my $id (@ids) {
	my $obj = $db->get_Seq_by_id($id);
  	my $seq = $obj->seq;
	my $chr_length  = $obj->length;
	my $chr_reads = int(($chr_length * $fraction)/$read_len);

	if ( $id =~ /^Y$/ ) {
		$chr_reads = 0;
	}

	my $SCE_start = 0;

	if ( grep(/^$id/, @SCE) ) {
		my $lower = 1;
		my $upper = $chr_length;
		$SCE_start = int($lower+rand($upper - $lower));
		
		if (exists($acrocent_chroms{$id})) {
			while ($SCE_start < $acrocent_chroms{$id}->[1]+1_000_000) {
				$SCE_start = int($lower+rand($upper - $lower));
			}
		}

		print OUT "$id $SCE_start $prefix\n";
	}

	
	$created_reads += $chr_reads;

	warn "Processing $id\t$chr_length\t$chr_reads ...\n";

	my $lower = 1;
	my $upper = $chr_length;

	my $back_reads = $chr_reads*$background;
	my @read_nums = 1..$chr_reads;
	my @shuffled_indexes = shuffle(0..$#read_nums);
	# Get just N of them.
	my @pick_indexes = @shuffled_indexes[ 0 .. $back_reads - 1 ];


	foreach my $read_num (1..$chr_reads) { 
		$counter++;

		my $frg_length = gauss($mean_frg_length, $SD_frg_length);

		my $randA = int(rand($chr_length-($frg_length-1))) + 1;
		redo if $seen{$randA}++;

		my $endA = $randA+$frg_length;

		my $subseqA = $obj->subseq($randA => $endA);
		my $pattern = 'N'x($read_len/2);		
		redo if $subseqA =~ /$pattern/;

		my $read;
		my $todo = rand();
		
		my $frw_start = $randA;
		my $rev_start = ($endA - $read_len);

		my $begin = substr($subseqA, 0, $read_len) if $todo < 0.5;
		my $end = substr($subseqA, -$read_len) if $todo > 0.5;

		if (defined($begin)) {
			$read = $begin;
		} else {
			$read = $end;
		}

		if (grep(/$read_num/, @pick_indexes)) {
			$read =~ tr/ACGTacgt/TGCAtgca/;
			$read = reverse($read);
		}

		if ( grep(/^$id/, @SCE) ) {

			if ( grep(/^$id/, @wc) ) {
				if ( ($frw_start >= $lower and $frw_start <= $SCE_start) or ($rev_start >= $lower and $rev_start <= $SCE_start) ) {
					$todo = rand();
					if ($todo > 0.5) {	
						$read =~ tr/ACGTacgt/TGCAtgca/;
						$read = reverse($read);
					}
				} else {		
					$read =~ tr/ACGTacgt/TGCAtgca/;
					$read = reverse($read);
					
				}
		
			} elsif ( grep(/$id/, @c) ) {
				if ( ($frw_start >= $lower and $frw_start <= $SCE_start) or ($rev_start >= $lower and $rev_start <= $SCE_start) ) {
					$read =~ tr/ACGTacgt/TGCAtgca/;
					$read = reverse($read);
				} else {
					$todo = rand();
					if ($todo > 0.5) {	
						$read =~ tr/ACGTacgt/TGCAtgca/;
						$read = reverse($read);
					}
				}
			} else {
				if ( ($frw_start >= $lower and $frw_start <= $SCE_start) or ($rev_start >= $lower and $rev_start <= $SCE_start) ) {
				} else {
					$todo = rand();
					if ($todo > 0.5) {	
						$read =~ tr/ACGTacgt/TGCAtgca/;
						$read = reverse($read);
					}
				}
			}
					

		} else {		
			if (grep(/^$id/, @wc)) {
				$todo = rand();
				if ($todo > 0.5) {	
					$read =~ tr/ACGTacgt/TGCAtgca/;
					$read = reverse($read);
				}
			} elsif (grep(/^$id/, @c)) {
				$read =~ tr/ACGTacgt/TGCAtgca/;
				$read = reverse($read);
			}
		}			

		print $read1 '@readA_'.$frg_length.'_'.$counter."\n";
		print $read1 "$read\n";
		print $read1 "+\n";
		print $read1 '2' x length($read),"\n";	
		
	}
}


## Align simulated reads using BWA
my $cmd;
warn "Aligning\n";
$cmd = join( ' ', 'bwa mem',
                  '-t', 4,
                  $GENOME,
		  './read1.fq',
                  '>', $prefix.'.sam'
               );
system $cmd;

warn "Sam2Bam\n";
$cmd = join( ' ', 'samtools', 'view',
                 '-O BAM',
                  $prefix.'.sam',
                  '-o', $prefix.'.bam',
               );
system $cmd;

warn "Sort\n";
$cmd = join( ' ', 'samtools', 'sort',
                  '-O BAM',
                  $prefix.'.bam',
                  '-o', $prefix.'_srt.bam',
               );
system $cmd;

warn "Marking duplicates\n";
$cmd = join( ' ', 'sambamba', 'markdup',
                  $prefix.'_srt.bam',
		  $prefix.'_srt_dedup.bam',
               );
system $cmd;

warn "Index\n";
$cmd = join( ' ', 'samtools', 'index',
                  $prefix.'_srt_dedup.bam',
               );
system $cmd;

warn "Clean\n";
$cmd = join( ' ', 'rm', $prefix.'.sam',
			$prefix.'.bam',
			$prefix.'_srt.bam',
		);
system $cmd;

srand();

close $read1;
}

## Subroutines ##
#################
## Simulate Gauss distribution of insert sizes
sub gauss {	
	# argumetns: mean insert size and SD
	my ($mean,$sd) = @_;
	my $rate1;
	my $rate2;
	my $random1;
	my $below = $mean-3*$sd;
	my $above = $mean+3*$sd;

		$rate1 = rand(1);
		$rate2 = rand(1);
		$random1 = sqrt(-2*log($rate1))*sin(2*PI*$rate2);
		$random1 = int ($mean+($random1*$sd));
		#($random1>$below) && ($random1<$above) && return $random1;
		return $random1;
}

## Format exported fasta sequence
sub to_fasta {
        my ($seqName, $seq, $len) = @_;

	# default to 80 characters of sequence per line
        $len = 80 unless $len;

        my $formatted_seq = ">$seqName\n";
        while (my $chunk = substr($seq, 0, $len, "")) {
                $formatted_seq .= "$chunk\n";
        }

        return $formatted_seq;
}


## Print help
sub help { 
print "\nsimulateStrandSreads.pl <reference genome> <fraction> <background> <index>\n\n";
print "<reference genome> : FASTA file of the reference genome to sample reads from.\n";
print "<fraction> : Fraction of the reference genome to sample reads from.\n";
print "<background> : Level of noise to simulate.\n";
print "<index> : Unique ID for a given set of simulated libraries.\n\n";
}

