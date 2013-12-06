#!/usr/bin/perl/ -w
$|++;
use strict;
use File::Path;
use Time::HiRes qw( time );

######################################################################################################################################################
#
#	Description
#		This is a perl script to label the multiple hits in a sam file (i.e. add the NH:i value and change non-unique hit to 0 mapping quality). 
#	Unaligned sequence will be output in a new file in the format of fastq.
#
#	Input	
#
#		--samPath=			the samPath
#		--maxHitQual=		the maximum number of hit that will be given a mapping quality of 255, anything larger than this number will have a quality of 0  
#		--maxHitRetain= 	the maximum number of hit that will be retained;
#		--sepPEUnaligned= 	"yes" or "no", to seperate the unaligned PE seqeunces into two mate-files; default = "no";
#		--outDir= 			path of the output dir
#		--rmOriginalSam= 	"yes" or "no". Remove the original sam file or not. default no.
#		--shortenRdNameBol=	"yes" or "no". shorten the read name or not;
#		--shortenRdNameTag=	string; if shortenRdNameBol = 'yes', read names will be changed as shortenRdNameTag plus a number;
#		--header=			"yes" or "no". To inheritate the header from the original sam file or not.
#
#	Output
#		1. a modified sam file with NH:i number and modified mapping quality value (mainly for display purposes in IGV, i.e. value 0 appears as transparent)
#		2. a fastq file with unmapped reads, mainly for mapping of junctions in tophat and denovo assembly
#
#	Usage
#		perl SAMMultiHitAndUnalignLabeller_v0.5.pl --samPath=EhistolyticaGenomic_AmoebaDB-1.0_all_pair.100nt.1M_frag.250.20.300.200.sam --maxHitQual=1 --maxHitRetain=50 --sepPEUnaligned=no --rmOriginalSam=no --header=yes
#
#
#	Assumption 
#		1. Reads of the same name are continous in the sam file, which is the case in bowtie output, i.e. unsorted
#
#	Version history
#
#		v0.2
#		-Add the unaligned read num count;	
#		-corrected totalReadNum as real total read numbe, which was treated actually aligned read number;
#
#		v0.3
#		-SAM flag translation was added instead of just using 0, 4 or 16
#		-counting of NH number is now different as pairend reads always have two reads of the same name;
#		-mismatch option disabled, due the the potential difficulties in handling pairend read;
#		-random function turned off due to the tectnical difficulties in handling pairend read;
#		-speed improved by avoid passing variables into a highly repetitive subroutine
#
#		v0.4
#		-will output the unaligned fastq according to single or pair end reads;
#
#		v0.5
#		-the output dir option added;
#		-remove original SAM file option added;
#
#		v0.6
#		-added shortenRdNameTag option;
#
######################################################################################################################################################
#
#		Example sam file head
#
#@HD	VN:1.0	SO:unsorted
#@SQ	SN:DS571145	LN:530629
#@PG	ID:Bowtie	VN:0.12.7	CL:"bowtie -p 6 -k 40 -m 40 -n 2 -f -S --mapq 1 contig.DS571145_1000000rd_50nt DS571145_1000000rd_50nt.fas"
#EHI_151410_0_r7448_259_299_DS571145_54231_54271_cnt_50M_+	4	*	0	0	*	*	0	0	TAATTAATAAAATACTTTAGTTCTTAAAGAGAGTCATGAATNNNNNNNNN	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	XM:i:0
#EHI_153780_0_r45520_6437_6486_DS571145_509278_509327_cnt_50M_-	16	DS571145	509278	1	50M	*	0	0	CTTTTTCCTCCAAACTCTTAATTAATAACACAAACATTAACATTGCTTTA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	XA:i:0	MD:Z:50	NM:i:0
#EHI_153130_0_r84591_201_250_DS571145_391161_391210_cnt_50M_+	0	DS571145	391161	1	50M	*	0	0	AAATTTGATTATTTAGCAATAAGAAATTTAATCATTACTTTAAGACAAAT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	XA:i:0	MD:Z:50	NM:i:0

#==========================================================Main body starts==========================================================================#


#1----------Read the parameters----------#
use vars qw ($samPath $maxHitQual $maxHitRetain $sepPEUnaligned $rmOriginalSam $header $shortenRdNameBol $shortenRdNameTag $outDir);
($samPath, $maxHitQual, $maxHitRetain, $sepPEUnaligned, $rmOriginalSam, $header, $shortenRdNameBol, $shortenRdNameTag, $outDir) = readParameters();
printCMDLogOrFinishMessage("CMDLog");

#2----------On-the-fly edit the sam file----------#
onTheFlyEditSam();

#3----------remove the originalSam if the user choose to----------#
if ($rmOriginalSam eq "yes") {
	system ("rm $samPath");
	my $samFileName = filePathToFileName($samPath);
	print $samFileName." removed.\n";
}

printCMDLogOrFinishMessage("finishMessage");

exit;
#========================================================= Main body ends ===========================================================================#

################################################################## readParameters ####################################################################
sub readParameters {

	$rmOriginalSam = "no";
	$sepPEUnaligned = "no"; #---for HMMSplicer
	$header = "yes";
	
	foreach (@ARGV) {
		if ($_ =~ m/--samPath=/) {$samPath = substr ($_, index ($_, "=")+1);} 
		elsif ($_ =~ m/--maxHitQual=/) {$maxHitQual = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--maxHitRetain=/) {$maxHitRetain = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--sepPEUnaligned=/) {$sepPEUnaligned = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--rmOriginalSam=/) {$rmOriginalSam = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--header=/) {$header = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--shortenRdNameBol=/) {$shortenRdNameBol = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--shortenRdNameTag=/) {$shortenRdNameTag = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--outDir=/) {$outDir = substr ($_, index ($_, "=")+1);}
	}
	
	system "mkdir -p -m 777 $outDir";
	
	return ($samPath, $maxHitQual, $maxHitRetain, $sepPEUnaligned, $rmOriginalSam, $header, $shortenRdNameBol, $shortenRdNameTag, $outDir);
}
#################################################### onTheFlyEditSam #################################################################################
sub onTheFlyEditSam {

	my $totalStart = time();

    print "Estimating the number of alignments in the sam file.\n";
	my (%filteredHitCountHsh, %originalHitCountHsh, %mismatchCountHsh);
	
	my @samPathSplt = split /\//, $samPath;
	my $outTag = $samPathSplt[-1];
	$outTag =~ s/\.sam//;

	#---estimate the number of lines in the file
	open (INFILE, $samPath) || die "Can't open $samPath.\n";
	my $tmpFile = $outTag."_tmp.txt";
	system "tail -100000 $samPath >$tmpFile";
	my $samPathSize = -s "$samPath";
   	my $tmpFileSize = -s "$tmpFile";
	system "rm $tmpFile";
   	my $samTotalLineNum = int (($samPathSize/$tmpFileSize)*100000);
   	print "Estimated to have ".$samTotalLineNum." alingments in the sam file.\n";

	#---read the sam
	open (EDITSAM, ">$outDir/$outTag.edited.sam");
	open (UNALGNFQSE, ">$outDir/$outTag.unaligned.fq");
	
	my $secondMateOut;
	if ($sepPEUnaligned eq "yes") {
		open (UNALGNFQPE1, ">$outDir/$outTag.unaligned.PE.1.fq");
		open (UNALGNFQPE2, ">$outDir/$outTag.unaligned.PE.2.fq");
		$secondMateOut = *UNALGNFQPE2;
	} else {
		open (UNALGNFQPE1, ">$outDir/$outTag.unaligned.fq");
		$secondMateOut = *UNALGNFQPE1;
	}

	#---get rid of the header and get the 1st read
	my ($theCurntLine, $theNextLine);
	while (my $theLine = <INFILE>) {
		chomp $theLine;
		if (not($theLine =~ m/^\@/)) {
			$theCurntLine = $theLine;
			last;
		} else {#---copy the header
			print EDITSAM $theLine."\n" if ($header eq "yes");
		}
	}
	
	#---go through all the lines
	my ($outputLine, $changeRead);
	my @outputAry;
	my $lineProc = 0;
	my $progCount = 0;
	my $totalReadNum = 0;
	my $originalAlignNum = 0;
	my $passedReadNum = 0;
	my $passedAlignNum = 0;
	my $mismatchFiltered = 0;
	my $hitFiltered = 0;
	my $tmpOriginalHitCount = 0;
	my %hitNumIndexHsh;
	my $unalignedReadNum = 0;
	my $alignedReadNum = 0;

	print "Start process the alignments.\n"; 
	my ($intervalStart, $intervalEnd);
	$intervalStart = time();
	
	#---define the SAM Table
	my $SAMFlagTableHsh_ref = defineSAMFlagTable();
	my %SAMFlagTableHsh = %{$SAMFlagTableHsh_ref};
	
	#--- define hashes to temporarily hold the unaligned pairend read
	my (%tmpPESeqHsh, %tmpPEQualHsh);
	
	while ($theNextLine = <INFILE>) {
		$lineProc++;
		$progCount++;
		
		if ($progCount == 100000) {
			print "$lineProc aligns processed. "; 
			$progCount=0;
			$intervalEnd = time();
			my $timeElapsed = $intervalEnd - $intervalStart;
			$timeElapsed = sprintf ("%.2f", $timeElapsed);
			my $estimatedEnd = (($samTotalLineNum - $lineProc)*$timeElapsed)/100000;
			$estimatedEnd = sprintf ("%.2f", $estimatedEnd/60);
			print "Time for last 100000 aligns is ".$timeElapsed." sec. Estimated to end in ".$estimatedEnd." mins.\n";
			$intervalStart = time();
		}
		
		chomp $theNextLine;
		my @theCurntLineSplt = split (/\t/, $theCurntLine);
		my @theNextLineSplt = split (/\t/, $theNextLine);
		
		my $theCurntLineSplt_ref = \@theCurntLineSplt;
		my $theNextLineSplt_ref = \@theNextLineSplt;
		
		#---store all SAM bits in a Hsh
		my $SAMFlag = $theCurntLineSplt[1];
		my $SAMBitStr = $SAMFlagTableHsh{$SAMFlag};
		my @SAMBitAry = split /\+/, $SAMBitStr;
		my %SAMBitHsh;
		foreach my $SAMBit (@SAMBitAry) {$SAMBitHsh{$SAMBit}++;}
		
		if (exists $SAMBitHsh{4}) {#---unaligned reads as the bit 4 is present

			$originalHitCountHsh{0}++;
			if (exists $SAMBitHsh{1}) {#---PEnd
				my $firstSecond;
				$firstSecond = 1 if $SAMBitHsh{64};
				$firstSecond = 2 if $SAMBitHsh{128};
				${$tmpPESeqHsh{$theCurntLineSplt[0]}}{$firstSecond} = $theCurntLineSplt[9];
				${$tmpPEQualHsh{$theCurntLineSplt[0]}}{$firstSecond} = $theCurntLineSplt[10];
				
				#---print if pairs are both exsists
				if ((exists ${$tmpPESeqHsh{$theCurntLineSplt[0]}}{1}) and (exists ${$tmpPESeqHsh{$theCurntLineSplt[0]}}{2})) {
					print UNALGNFQPE1 "@".$theCurntLineSplt[0]."/1\n";
					print UNALGNFQPE1 ${$tmpPESeqHsh{$theCurntLineSplt[0]}}{1}."\n";
					print UNALGNFQPE1 "+".$theCurntLineSplt[0]."/1\n";
					print UNALGNFQPE1 ${$tmpPEQualHsh{$theCurntLineSplt[0]}}{1}."\n";
					print $secondMateOut "@".$theCurntLineSplt[0]."/2\n";
					print $secondMateOut ${$tmpPESeqHsh{$theCurntLineSplt[0]}}{2}."\n";
					print $secondMateOut "+".$theCurntLineSplt[0]."/2\n";
					print $secondMateOut ${$tmpPEQualHsh{$theCurntLineSplt[0]}}{2}."\n";
					delete $tmpPESeqHsh{$theCurntLineSplt[0]};
					delete $tmpPEQualHsh{$theCurntLineSplt[0]};
				}
				
			} else {#---SEnd
				print UNALGNFQSE "@".$theCurntLineSplt[0]."\n";
				print UNALGNFQSE $theCurntLineSplt[9]."\n";
				print UNALGNFQSE "+".$theCurntLineSplt[0]."\n";
				my $qual = $theCurntLineSplt[10];
				$qual =~ s/J/I/g;
				print UNALGNFQSE $qual."\n";
			}
			$totalReadNum++;
			$unalignedReadNum++;

		} else {#---aligned read as the bit 4 is absent

			push @outputAry, $theCurntLine; #---aligned, to be output

			$originalAlignNum++;
			$tmpOriginalHitCount++;
			
			#----execute the following commands no matter changeRead is yes or no
			if ($theCurntLineSplt[0] eq $theNextLineSplt[0]) {#---read ID not the same, i.e. unique read
				$changeRead = "no";
			} else {
				$changeRead = "yes";
				$totalReadNum++;
				$alignedReadNum++;
				$originalHitCountHsh{$tmpOriginalHitCount}++;
				$hitNumIndexHsh{$tmpOriginalHitCount} = "yes";
				$tmpOriginalHitCount = 0;
			}
			
			if ($changeRead eq "yes") {
				
				my $hitNum  = @outputAry; #---single end
				$hitNum = int ($hitNum/2) if (exists $SAMBitHsh{1});#---pair end if SAM Flag has 1
				
				if ($hitNum <= $maxHitRetain) {#---if the read passed the hitnum filter
					$passedReadNum++;
					$passedAlignNum += $hitNum;
					$hitNumIndexHsh{$hitNum} = "yes";
	
					my $qual = 255;
					$qual = 0 if ($hitNum > $maxHitQual);
					$filteredHitCountHsh{$hitNum}++;
					foreach my $alignmentToPrint (@outputAry) {
						my @SAMColumnAry = split /\t/, $alignmentToPrint;
						$SAMColumnAry[4] = $qual;
						my $NHStr = "NH:i:".$hitNum;
						push @SAMColumnAry, $NHStr;
						$SAMColumnAry[0] = $shortenRdNameTag.$passedReadNum if $shortenRdNameBol eq 'yes';
						my $outStr = join "\t", @SAMColumnAry;
						print EDITSAM $outStr."\n";
					}
					
				} else {
					$hitFiltered++;
				}

				@outputAry = (); #---empty the Ary
			}
		}

		if (eof(INFILE)) {#---the last line

			#---store all SAM bits in a Hsh
			my $SAMFlag = $theNextLineSplt[1];
			my $SAMBitStr = $SAMFlagTableHsh{$SAMFlag};
			my @SAMBitAry = split /\+/, $SAMBitStr;
			my %SAMBitHsh;
			foreach my $SAMBit (@SAMBitAry) {$SAMBitHsh{$SAMBit}++;}
			
			if (exists $SAMBitHsh{4}) {#---unaligned reads as the bit 4 is present

				$originalHitCountHsh{0}++;
	
				if (exists $SAMBitHsh{1}) {#---PEnd
					my $firstSecond;
					$firstSecond = 1 if $SAMBitHsh{64};
					$firstSecond = 2 if $SAMBitHsh{128};
					${$tmpPESeqHsh{$theNextLineSplt[0]}}{$firstSecond} = $theNextLineSplt[9];
					${$tmpPEQualHsh{$theNextLineSplt[0]}}{$firstSecond} = $theNextLineSplt[10];
					
					#---print if pairs are both exsists
					if ((exists ${$tmpPESeqHsh{$theNextLineSplt[0]}}{1}) and (exists ${$tmpPESeqHsh{$theNextLineSplt[0]}}{2})) {
						print UNALGNFQPE1 "@".$theNextLineSplt[0]."/1\n";
						print UNALGNFQPE1 ${$tmpPESeqHsh{$theNextLineSplt[0]}}{1}."\n";
						print UNALGNFQPE1 "+".$theNextLineSplt[0]."/1\n";
						print UNALGNFQPE1 ${$tmpPEQualHsh{$theNextLineSplt[0]}}{1}."\n";
						print $secondMateOut "@".$theNextLineSplt[0]."/2\n";
						print $secondMateOut ${$tmpPESeqHsh{$theNextLineSplt[0]}}{2}."\n";
						print $secondMateOut "+".$theNextLineSplt[0]."/2\n";
						print $secondMateOut ${$tmpPEQualHsh{$theNextLineSplt[0]}}{2}."\n";
						delete $tmpPESeqHsh{$theNextLineSplt[0]};
						delete $tmpPEQualHsh{$theNextLineSplt[0]};
					}
					
				} else {#---SEnd
				
					print UNALGNFQSE "@".$theNextLineSplt[0]."\n";
					print UNALGNFQSE $theNextLineSplt[9]."\n";
					print UNALGNFQSE "+".$theNextLineSplt[0]."\n";
					print UNALGNFQSE $theNextLineSplt[10]."\n";
				}
	
				$totalReadNum++;
				$unalignedReadNum++;

			} else {
				
				push @outputAry, $theNextLine;

				$originalAlignNum++;
				$tmpOriginalHitCount++;

				$totalReadNum++;
				$alignedReadNum++;
				$originalHitCountHsh{$tmpOriginalHitCount}++;
				$hitNumIndexHsh{$tmpOriginalHitCount} = "yes";
				$tmpOriginalHitCount = 0;
			
				my $hitNum  = @outputAry; #---single end
				$hitNum = int ($hitNum/2) if (exists $SAMBitHsh{1});#---pair end if SAM Flag has 1
				
				if ($hitNum <= $maxHitRetain) {#---if the read passed the hitnum filter
					$passedReadNum++;
					$passedAlignNum += $hitNum;
					$hitNumIndexHsh{$hitNum} = "yes";
	
					my $qual = 255;
					$qual = 0 if ($hitNum > $maxHitQual);
					$filteredHitCountHsh{$hitNum}++;
					foreach my $alignmentToPrint (@outputAry) {
						my @SAMColumnAry = split /\t/, $alignmentToPrint;
						$SAMColumnAry[4] = $qual;
						my $NHStr = "NH:i:".$hitNum;
						push @SAMColumnAry, $NHStr;
						$SAMColumnAry[0] = $shortenRdNameTag.$passedReadNum if $shortenRdNameBol eq 'yes';
						my $outStr = join "\t", @SAMColumnAry;
						print EDITSAM $outStr."\n";
					}
				} else {
					$hitFiltered++;
				}
			}
		}

		#---next line becomes current line
		$theCurntLine = $theNextLine;

	}
	
	print "$lineProc alignments processed.\n"; 
	my $totalEnd = time();
	my $timeElapsedSec = $totalEnd-$totalStart;
	my $timeElapsedSecPer100000Align = sprintf ("%.2f", ($timeElapsedSec/$originalAlignNum)*100000);
	my $timeElapsedMin = sprintf ("%.2f", $timeElapsedSec/60);
	
	print "\nTotal time elapsed = ".$timeElapsedMin." minutes\n";
	print "Average time elapsed per 100000 alignments = ".$timeElapsedSecPer100000Align." seconds\n";

	close INFILE;
	close EDITSAM;
	close UNALGNFQSE;

	#---print the statistics
	open (LOGFILE, ">$outDir/$outTag.stats.log");
	print LOGFILE "Total time tlapsed = ".$timeElapsedMin." minutes\n";
	print LOGFILE "Average time tlapsed per 100000 alignments = ".$timeElapsedSecPer100000Align." seconds\n\n";
	print LOGFILE "totalReadNum = ".$totalReadNum."\n";
	my $unalignedReadPct = sprintf ("%.2f", ($unalignedReadNum/$totalReadNum)*100);
	my $alignedReadPct = sprintf ("%.2f", ($alignedReadNum/$totalReadNum)*100);
	print LOGFILE "unalignedReadNum = ".$unalignedReadNum."\t(".$unalignedReadPct."%)\n";
	print LOGFILE "alignedReadNum = ".$alignedReadNum."\t(".$alignedReadPct."%)\n";
	print LOGFILE "passedReadNum = ".$passedReadNum."\n";
	print LOGFILE "originalAlignNum = ".$originalAlignNum."\n";
	print LOGFILE "passedAlignNum = ".$passedAlignNum."\n";
	print LOGFILE "number of alignments filtered by mismatch = ".$mismatchFiltered."\n";
	print LOGFILE "number of read filtered by hitNum= ".$hitFiltered."\n";
	print LOGFILE "\nmismatch\tAlignmentNum\tAlignmentPct\n";
	foreach my $mismatch (sort {$a <=> $b} keys %mismatchCountHsh) {
		my $alignmentNum = $mismatchCountHsh{$mismatch};
		my $alignmentPct = sprintf ("%.4f", ($alignmentNum/$originalAlignNum)*100);
		print LOGFILE $mismatch."\t".$alignmentNum."\t".$alignmentPct."\n";
	}

	print LOGFILE "\nhitNum\toriginalNum\tfilteredNum\toriginalPct\tfilteredPct\n";
	foreach my $hitNum (sort {$a <=> $b} keys %hitNumIndexHsh) {
		my ($filteredNum, $filteredPct, $originalNum, $originalPct);
		
		if (exists $filteredHitCountHsh{$hitNum}) {
			$filteredNum = $filteredHitCountHsh{$hitNum};
			$filteredPct = sprintf ("%.4f", ($filteredNum/$passedReadNum)*100);
		} else {
			$filteredNum = $filteredPct = "null";
		}
		
		if (exists $originalHitCountHsh{$hitNum}) {
			$originalNum = $originalHitCountHsh{$hitNum};
			$originalPct = sprintf ("%.4f", ($originalNum/$totalReadNum)*100);
		} else {
			$originalNum = $originalPct = "null";
		}
		
		print LOGFILE $hitNum."\t".$originalNum."\t".$filteredNum."\t".$originalPct."%\t".$filteredPct."%\n";
	}
	close LOGFILE;
}
##################################################################################################################################################
sub defineSAMFlagTable {
	
#
#copied from http://bioinformatics.bc.edu/chuanglab/wiki/index.php/SAM_pairwise_flag_translation_table
#
#0x0001 1 the read is paired in sequencing, no matter whether it is mapped in a pair 
#0x0002 2 the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment)  
#0x0004 4 the query sequence itself is unmapped 
#0x0008 8 the mate is unmapped  
#0x0010 16 strand of the query (0 for forward; 1 for reverse strand) 
#0x0020 32 strand of the mate  
#0x0040 64 the read is the ﬁrst read in a pair  
#0x0080 128 the read is the second read in a pair 
#0x0100 256 the alignment is not primary (a read having split hits may have multiple primary alignment records) 
	
	my %SAMFlagTableHsh = (

		0 => "0",
		1 => "1",
		2 => "2",
		3 => "1+2",
		4 => "0+4",
		5 => "1+4",
		6 => "0+2+4",
		7 => "1+2+4",
		8 => "0+8",
		9 => "1+8",
		10 => "0+2+8",
		11 => "1+2+8",
		12 => "0+4+8",
		13 => "1+4+8",
		14 => "0+2+4+8",
		15 => "1+2+4+8",
		16 => "0+16",
		17 => "1+16",
		18 => "0+2+16",
		19 => "1+2+16",
		20 => "0+4+16",
		21 => "1+4+16",
		22 => "0+2+4+16",
		23 => "1+2+4+16",
		24 => "0+8+16",
		25 => "1+8+16",
		26 => "0+2+8+16",
		27 => "1+2+8+16",
		28 => "0+4+8+16",
		29 => "1+4+8+16",
		30 => "0+2+4+8+16",
		31 => "1+2+4+8+16",
		32 => "0+32",
		33 => "1+32",
		34 => "0+2+32",
		35 => "1+2+32",
		36 => "0+4+32",
		37 => "1+4+32",
		38 => "0+2+4+32",
		39 => "1+2+4+32",
		40 => "0+8+32",
		41 => "1+8+32",
		42 => "0+2+8+32",
		43 => "1+2+8+32",
		44 => "0+4+8+32",
		45 => "1+4+8+32",
		46 => "0+2+4+8+32",
		47 => "1+2+4+8+32",
		48 => "0+16+32",
		49 => "1+16+32",
		50 => "0+2+16+32",
		51 => "1+2+16+32",
		52 => "0+4+16+32",
		53 => "1+4+16+32",
		54 => "0+2+4+16+32",
		55 => "1+2+4+16+32",
		56 => "0+8+16+32",
		57 => "1+8+16+32",
		58 => "0+2+8+16+32",
		59 => "1+2+8+16+32",
		60 => "0+4+8+16+32",
		61 => "1+4+8+16+32",
		62 => "0+2+4+8+16+32",
		63 => "1+2+4+8+16+32",
		64 => "0+64",
		65 => "1+64",
		66 => "0+2+64",
		67 => "1+2+64",
		68 => "0+4+64",
		69 => "1+4+64",
		70 => "0+2+4+64",
		71 => "1+2+4+64",
		72 => "0+8+64",
		73 => "1+8+64",
		74 => "0+2+8+64",
		75 => "1+2+8+64",
		76 => "0+4+8+64",
		77 => "1+4+8+64",
		78 => "0+2+4+8+64",
		79 => "1+2+4+8+64",
		80 => "0+16+64",
		81 => "1+16+64",
		82 => "0+2+16+64",
		83 => "1+2+16+64",
		84 => "0+4+16+64",
		85 => "1+4+16+64",
		86 => "0+2+4+16+64",
		87 => "1+2+4+16+64",
		88 => "0+8+16+64",
		89 => "1+8+16+64",
		90 => "0+2+8+16+64",
		91 => "1+2+8+16+64",
		92 => "0+4+8+16+64",
		93 => "1+4+8+16+64",
		94 => "0+2+4+8+16+64",
		95 => "1+2+4+8+16+64",
		96 => "0+32+64",
		97 => "1+32+64",
		98 => "0+2+32+64",
		99 => "1+2+32+64",
		100 => "0+4+32+64",
		101 => "1+4+32+64",
		102 => "0+2+4+32+64",
		103 => "1+2+4+32+64",
		104 => "0+8+32+64",
		105 => "1+8+32+64",
		106 => "0+2+8+32+64",
		107 => "1+2+8+32+64",
		108 => "0+4+8+32+64",
		109 => "1+4+8+32+64",
		110 => "0+2+4+8+32+64",
		111 => "1+2+4+8+32+64",
		112 => "0+16+32+64",
		113 => "1+16+32+64",
		114 => "0+2+16+32+64",
		115 => "1+2+16+32+64",
		116 => "0+4+16+32+64",
		117 => "1+4+16+32+64",
		118 => "0+2+4+16+32+64",
		119 => "1+2+4+16+32+64",
		120 => "0+8+16+32+64",
		121 => "1+8+16+32+64",
		122 => "0+2+8+16+32+64",
		123 => "1+2+8+16+32+64",
		124 => "0+4+8+16+32+64",
		125 => "1+4+8+16+32+64",
		126 => "0+2+4+8+16+32+64",
		127 => "1+2+4+8+16+32+64",
		128 => "0+128",
		129 => "1+128",
		130 => "0+2+128",
		131 => "1+2+128",
		132 => "0+4+128",
		133 => "1+4+128",
		134 => "0+2+4+128",
		135 => "1+2+4+128",
		136 => "0+8+128",
		137 => "1+8+128",
		138 => "0+2+8+128",
		139 => "1+2+8+128",
		140 => "0+4+8+128",
		141 => "1+4+8+128",
		142 => "0+2+4+8+128",
		143 => "1+2+4+8+128",
		144 => "0+16+128",
		145 => "1+16+128",
		146 => "0+2+16+128",
		147 => "1+2+16+128",
		148 => "0+4+16+128",
		149 => "1+4+16+128",
		150 => "0+2+4+16+128",
		151 => "1+2+4+16+128",
		152 => "0+8+16+128",
		153 => "1+8+16+128",
		154 => "0+2+8+16+128",
		155 => "1+2+8+16+128",
		156 => "0+4+8+16+128",
		157 => "1+4+8+16+128",
		158 => "0+2+4+8+16+128",
		159 => "1+2+4+8+16+128",
		160 => "0+32+128",
		161 => "1+32+128",
		162 => "0+2+32+128",
		163 => "1+2+32+128",
		164 => "0+4+32+128",
		165 => "1+4+32+128",
		166 => "0+2+4+32+128",
		167 => "1+2+4+32+128",
		168 => "0+8+32+128",
		169 => "1+8+32+128",
		170 => "0+2+8+32+128",
		171 => "1+2+8+32+128",
		172 => "0+4+8+32+128",
		173 => "1+4+8+32+128",
		174 => "0+2+4+8+32+128",
		175 => "1+2+4+8+32+128",
		176 => "0+16+32+128",
		177 => "1+16+32+128",
		178 => "0+2+16+32+128",
		179 => "1+2+16+32+128",
		180 => "0+4+16+32+128",
		181 => "1+4+16+32+128",
		182 => "0+2+4+16+32+128",
		183 => "1+2+4+16+32+128",
		184 => "0+8+16+32+128",
		185 => "1+8+16+32+128",
		186 => "0+2+8+16+32+128",
		187 => "1+2+8+16+32+128",
		188 => "0+4+8+16+32+128",
		189 => "1+4+8+16+32+128",
		190 => "0+2+4+8+16+32+128",
		191 => "1+2+4+8+16+32+128",
		192 => "0+64+128",
		193 => "1+64+128",
		194 => "0+2+64+128",
		195 => "1+2+64+128",
		196 => "0+4+64+128",
		197 => "1+4+64+128",
		198 => "0+2+4+64+128",
		199 => "1+2+4+64+128",
		200 => "0+8+64+128",
		201 => "1+8+64+128",
		202 => "0+2+8+64+128",
		203 => "1+2+8+64+128",
		204 => "0+4+8+64+128",
		205 => "1+4+8+64+128",
		206 => "0+2+4+8+64+128",
		207 => "1+2+4+8+64+128",
		208 => "0+16+64+128",
		209 => "1+16+64+128",
		210 => "0+2+16+64+128",
		211 => "1+2+16+64+128",
		212 => "0+4+16+64+128",
		213 => "1+4+16+64+128",
		214 => "0+2+4+16+64+128",
		215 => "1+2+4+16+64+128",
		216 => "0+8+16+64+128",
		217 => "1+8+16+64+128",
		218 => "0+2+8+16+64+128",
		219 => "1+2+8+16+64+128",
		220 => "0+4+8+16+64+128",
		221 => "1+4+8+16+64+128",
		222 => "0+2+4+8+16+64+128",
		223 => "1+2+4+8+16+64+128",
		224 => "0+32+64+128",
		225 => "1+32+64+128",
		226 => "0+2+32+64+128",
		227 => "1+2+32+64+128",
		228 => "0+4+32+64+128",
		229 => "1+4+32+64+128",
		230 => "0+2+4+32+64+128",
		231 => "1+2+4+32+64+128",
		232 => "0+8+32+64+128",
		233 => "1+8+32+64+128",
		234 => "0+2+8+32+64+128",
		235 => "1+2+8+32+64+128",
		236 => "0+4+8+32+64+128",
		237 => "1+4+8+32+64+128",
		238 => "0+2+4+8+32+64+128",
		239 => "1+2+4+8+32+64+128",
		240 => "0+16+32+64+128",
		241 => "1+16+32+64+128",
		242 => "0+2+16+32+64+128",
		243 => "1+2+16+32+64+128",
		244 => "0+4+16+32+64+128",
		245 => "1+4+16+32+64+128",
		246 => "0+2+4+16+32+64+128",
		247 => "1+2+4+16+32+64+128",
		248 => "0+8+16+32+64+128",
		249 => "1+8+16+32+64+128",
		250 => "0+2+8+16+32+64+128",
		251 => "1+2+8+16+32+64+128",
		252 => "0+4+8+16+32+64+128",
		253 => "1+4+8+16+32+64+128",
		254 => "0+2+4+8+16+32+64+128",
		255 => "1+2+4+8+16+32+64+128",
	);
	 
	 return \%SAMFlagTableHsh;
}
#################################################### print command log ##################################################################################
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptNameXext = $0;
		$scriptNameXext =~ s/\.\w+$//;
		open (CMDLOG, ">>$scriptNameXext.cmd.log.txt"); #---append the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}
	
}
#################################################### formatedDateTime ##################################################################################
sub filePathToFileName {

	my $filePath = $_[0];
	my @filePathSplt = split /\//, $filePath;
	my $fileName = $filePathSplt[-1];
	
	return $fileName;
}
