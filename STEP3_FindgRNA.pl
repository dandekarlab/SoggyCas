#!/usr/bin/perl -w
use strict;
use DBI;
use Bio::Restriction::EnzymeCollection;

my $resultfile1 = "TableData\/gRNATable.txt";
my $resultfile2 = "TableData\/20MerTable.txt";
my $resultfile3 = "TableData\/Temp_12MerTable.txt";
my $resultfile4 = "TableData\/Temp_8MerTable.txt";

my $conn = DBI->connect("dbi:mysql:soggycasdb",'root','') or die "Can't establish DB connectivity\n";

my %RestrictionEnzyme = ();
my $all_collection = Bio::Restriction::EnzymeCollection->new();
my $six_cutter_collection1 = $all_collection->cutters(4);
my $six_cutter_collection2 = $all_collection->cutters(6);
for my $enz( $six_cutter_collection1->each_enzyme())
{
	my $enzyme = $enz->name;
	my $site = $enz->site;
	$site=~tr/^//d;
	if($site=~/[^ATGC]/gi)
	{
		next;
	}
	$RestrictionEnzyme{$site}=$enzyme;
}

for my $enz( $six_cutter_collection2->each_enzyme())
{
	my $enzyme = $enz->name;
	my $site = $enz->site;
	$site=~tr/^//d;
	if($site=~/[^ATGC]/gi)
	{
		next;
	}
	$RestrictionEnzyme{$site} = $enzyme;
}

open OUT1,">$resultfile1" or die "can't open $resultfile1 for writing\n";
open OUT2,">$resultfile2" or die "can't open $resultfile2 for writing\n";
open OUT3,">$resultfile3" or die "can't open $resultfile3 for writing\n";
open OUT4,">$resultfile4" or die "can't open $resultfile4 for writing\n";
my $query = "select geneid,geneseq,genestart,chromosome from gene where geneid= 108996675";
my $sth = $conn->prepare($query);
$sth->execute();
my $geneid = '';
my $geneseq = '';
my $genestart = '';
my $chromosome = '';
my %gRNAIds = ();
my $GRNACnt = 0;
my %TwlMerInfo = ();
my %EghtMerInfo = ();
$sth->bind_columns(undef,\$geneid,\$geneseq,\$genestart,\$chromosome);
while($sth->fetch())
{
	&findgRNA($geneseq,$geneid,$chromosome,$genestart,\%gRNAIds,\$GRNACnt,\%TwlMerInfo,\%EghtMerInfo,\*OUT1,\*OUT2);
}
$sth->finish();
$conn->disconnect();
close OUT1;
close OUT2;

foreach my $eachTwlMer(keys %TwlMerInfo)
{
	my $gids = join",",@{$TwlMerInfo{$eachTwlMer}};
	print OUT3 "$eachTwlMer\t$gids\n";
}
close OUT3;

foreach my $eachEgtMer(keys %EghtMerInfo)
{
	my $gids = join",",@{$EghtMerInfo{$eachEgtMer}};
	print OUT4 "$eachEgtMer\t$gids\n";
}

close OUT4;


sub findgRNA
{
	my $seq = shift;
	my $gid = shift;
	my $chr = shift;
	my $gSt = shift;
	my $ref_gRNAIds = shift;
	my $ref_GRNACnt = shift;
	my $RefTwlMerInfo = shift;
	my $RefEghtMerInfo = shift;
	my $refOUT1 = shift;
	my $refOUT2 = shift;
	
	$seq=~tr/ //d;
	$seq=~s/\r//g;
	my $seqLen = length($seq);
	my %CheckMer = ();
	#print "GRNA\tStrand\tGC Percentage\tPosition\tMelTemp\t20+PAM N.Targets\t12+PAM N.Targets\t8+PAM N.Targets\tR Enzyme(s)\n";
	foreach  my $pos(0..($seqLen-1))
	{
		my $gRNA = substr($seq,$pos,23);
		if(length($gRNA)==23)
		{
			my $tttFlag = 0;
			if($gRNA=~/.*GG$/ig)
			{
				my $subtar =substr($gRNA,0,20);
				my $gc = $subtar=~tr/GC//;
				my $gcPer = ($gc*100/20);
				my $stPos = $pos+1;
				my $edPos = $pos+20;
				my $melTem = 64.9 +(41*(($gc-16.4)/20));
				my $TwlMer = substr($gRNA,8,15);
				my $EightMer = substr($gRNA,12,12);
				my $rEnzyme = '';
				foreach my $sites(keys %RestrictionEnzyme)
				{
					if($gRNA=~/$sites/)
					{
						if($rEnzyme eq '')
						{
							$rEnzyme=$RestrictionEnzyme{$sites};
						}
						else
						{
							$rEnzyme.=','.$RestrictionEnzyme{$sites};
						}
					}
				}
				my $guideRNAID = 'NA';
				my $gChrSt = $gSt+($stPos-1);
				my $gChrEd = $gSt+($edPos-1);
				$tttFlag =$gRNA=~s/TTTT/TTTT/g;
				$tttFlag = 0 if($tttFlag=~/^\s*$/);
				if(defined($$ref_gRNAIds{$gRNA}{'+'}))
				{
					$guideRNAID = $$ref_gRNAIds{$gRNA}{'+'};
					print $refOUT2 "$guideRNAID\t$gid\t$chr\t$gChrSt\t$gChrEd\n";
				}
				else
				{
					$$ref_GRNACnt++;
					$guideRNAID = 'GR'.$$ref_GRNACnt;
					$$ref_gRNAIds{$gRNA}{'+'} = $guideRNAID;
					print $refOUT1 "$guideRNAID\t$gRNA\t$gcPer\t$melTem\t$tttFlag\t$rEnzyme\n";
					print $refOUT2 "$guideRNAID\t$gid\t$chr\t$gChrSt\t$gChrEd\t\+\n";
				}
				if(!defined($CheckMer{$TwlMer}{$guideRNAID}))
				{
					push(@{$$RefTwlMerInfo{$TwlMer}},$guideRNAID);
					$CheckMer{$TwlMer}{$guideRNAID} = 0;
				}
				if(!defined($CheckMer{$EightMer}{$guideRNAID}))
				{
					push(@{$$RefEghtMerInfo{$EightMer}},$guideRNAID);
					$CheckMer{$EightMer}{$guideRNAID}=0;
				}
				
			}
			if($gRNA=~/^CC.*/ig)
			{
				my $subtar =substr($gRNA,3,20);
				my $gc = $subtar=~tr/GC//;
				my $gcPer = ($gc*100/20);
				my $stPos = $pos+4;
				my $edPos = $stPos+19;

				my $melTem = 64.9 +(41*(($gc-16.4)/20));
				my $TwlMer = substr($gRNA,0,15);
				my $EightMer = substr($gRNA,0,11);
				my $rEnzyme = '';
				foreach my $sites(keys %RestrictionEnzyme)
				{
					if($gRNA=~/$sites/)
					{
						if($rEnzyme eq '')
						{
							$rEnzyme=$RestrictionEnzyme{$sites};
						}
						else
						{
							$rEnzyme.=','.$RestrictionEnzyme{$sites};
						}
					}
				}
				my $guideRNAID = 'NA';
				my $gChrSt = $gSt+($stPos-1);
				my $gChrEd = $gSt+($edPos-1);
				$tttFlag =$gRNA=~s/TTTT/TTTT/g;
				$tttFlag = 0 if($tttFlag=~/^\s*$/);
				if(defined($$ref_gRNAIds{$gRNA}{'-'}))
				{
					$guideRNAID = $$ref_gRNAIds{$gRNA}{'-'};
					print $refOUT2 "$guideRNAID\t$gid\t$chr\t$gChrSt\t$gChrEd\n";
				}
				else
				{
					$$ref_GRNACnt++;
					$guideRNAID = 'GR'.$$ref_GRNACnt;
					$$ref_gRNAIds{$gRNA}{'-'} = $guideRNAID;
					print $refOUT1 "$guideRNAID\t$gRNA\t$gcPer\t$melTem\t$tttFlag\t$rEnzyme\n";
					print $refOUT2 "$guideRNAID\t$gid\t$chr\t$gChrSt\t$gChrEd\t\-\n";
				}
				if(!defined($CheckMer{$TwlMer}{$guideRNAID}))
				{
					push(@{$$RefTwlMerInfo{$TwlMer}},$guideRNAID);
					$CheckMer{$TwlMer}{$guideRNAID} = 0;
				}
				if(!defined($CheckMer{$EightMer}{$guideRNAID}))
				{
					push(@{$$RefEghtMerInfo{$EightMer}},$guideRNAID);
					$CheckMer{$EightMer}{$guideRNAID}=0;
				}
				
			}
		}
	}
}