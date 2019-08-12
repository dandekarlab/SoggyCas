#!/usr/bin/perl -w
use strict;
use Getopt::Long;

## Input Files ##
my $infile1 = ""; # GFF File #
my $infile2 = ""; # Protein Table #
my $infile3 = ""; # Genomic fasta file #

GetOptions ("gff=s" => \$infile1,
              "p=s"   => \$infile2,
              "fa=s"   => \$infile3
              )
or die("Error in command line arguments\n");


my $taxid = 51240;
my $organism = 'Juglans regia';
my $commonName = 'English walnut';

# Initialization #
my %mrnaInfo = ();
my %exonNo = ();
my %exonEnds = ();
my %intronNo = ();
my %mRNAMap = ();
my %ProteinInfo = ();
my %mrnaNos = ();
my %ProteinDetails = ();
my $TotalGenes = 0;
my $TotalmRNA = 0;
my $TotalProteins = 0;
my $totalExons = 0;
my $totalIntrons = 0;
my $TotalExonLen = 0;
my $TotalIntronLen = 0;

## Result files ##
my $resultfolder = 'TableData';
unless(-d $resultfolder)
{
	mkdir $resultfolder;
}
my $resulfile1 = "$resultfolder\/organismTable.txt";
my $resulfile2 = "$resultfolder\/geneTable.txt";
my $resulfile3 = "$resultfolder\/mrnaTable.txt";
my $resulfile4 = "$resultfolder\/exonTable.txt";
my $resulfile5 = "$resultfolder\/intronTable.txt";
my $resulfile6 = "$resultfolder\/proteinTable.txt";

my $header = '';
my $seq = '';
my %GenomeSeq = ();
open IN3,"<$infile3" or die "Can't open $infile3 for reading\n";
while(<IN3>)
{
	chomp;
	next if(/^\s*$/);
	next if(/^\s*\#+.*/);
	if(/^\s*\>+(\S+)\s*.*/)
	{
		if($header ne '')
		{
			$GenomeSeq{$header} = $seq;
			$header = '';
			$seq = '';
		}
		$header = $1;
	}
	else
	{
		$seq.=$_;
	}
}
close IN3;
if($header ne '')
{
	$GenomeSeq{$header} = $seq;
	$header = '';
	$seq = '';
}

## Reading GFF file and parsing contents ##
open IN1,"<$infile1" or die "Can't open $infile1 for reading\n";
open OUT2,">$resulfile2" or die "Can't open $resulfile2 for writing\n";
open OUT4,">$resulfile4" or die "Can't open $resulfile4 for writing\n";
open OUT5,">$resulfile5" or die "Can't open $resulfile5 for writing\n";
while(<IN1>)
{
	chomp;
	next if(/^\s*$/);
	next if(/^\s*\#+.*/);
	my @data = split"\t",$_,100;
	
	## Extracting Gene Information ##
	if($data[2] eq 'gene') 
	{
		my $chromosome = $data[0];
		my $geneStart = $data[3];
		my $geneEnd = $data[4];
		my $geneLength = ($geneEnd-$geneStart)+1;
		my $strand = $data[6];
		my $geneName = 'NA';
		my $bioType = 'NA';
		my $geneId = 'NA';
		my @info = split"\;",$data[8],100;
		foreach my $info1(@info)
		{
			my($ky,$vl) = split'\=',$info1;
			if($ky eq 'Name')
			{
				$geneName = $vl;
			}
			elsif($ky eq 'gene_biotype')
			{
				$vl=~tr/_/ /;
				$bioType = $vl;
			}
			elsif($ky eq 'Dbxref')
			{
				my @info1 = split"\:",$vl,100;
				if($info1[0] eq 'GeneID')
				{
					$geneId = $info1[$#info1];
				}
			}
		}
		$TotalGenes++;
		my $geneSeq = uc(substr($GenomeSeq{$chromosome},($geneStart-1),$geneLength));
		
		if($strand eq '-')
		{
			$geneSeq = reverse($geneSeq);
			$geneSeq=~tr/ATGC/TACG/;
		}
		print OUT2 "$taxid\t$chromosome\t$geneId\t$geneName\t$geneStart\t$geneEnd\t$geneLength\t$bioType\t$strand\t$geneSeq\n";
	}
	elsif($data[2] eq 'mRNA')
	{
		my $mrnaStart = $data[3];
		my $mrnaEnd = $data[4];
		my $mrnaId = 'NA';
		my $geneId = 'NA';
		my $Id = 'NA';
		my $mrnaLength = ($mrnaEnd-$mrnaStart)+1;
		my @info = split"\;",$data[8],100;
		foreach my $info1(@info)
		{
			my($ky,$vl) = split'\=',$info1;
			if($ky eq 'transcript_id')
			{
				$mrnaId = $vl;
			}
			if($ky eq 'ID')
			{
				$Id = $vl;
			}
			elsif($ky eq 'Dbxref')
			{
				my @info1 = split"\:",$vl,100;
				if($info1[0] eq 'GeneID')
				{
					$geneId = (split'\,',$info1[1])[0];
				}
			}
		}
		$mRNAMap{$Id} = $mrnaId;
		$mrnaNos{$mrnaId}++;
		$mrnaInfo{$mrnaNos{$mrnaId}}{$mrnaId} = "$geneId\t$mrnaId\t$mrnaStart\t$mrnaEnd\t$mrnaLength";
	}
	elsif($data[2] eq 'exon')
	{
		my $exonStart = $data[3];
		my $exonEnd = $data[4];
		my $exonLength = ($exonEnd-$exonStart)+1;
		my $mrnaId = 'NA';
		my @info = split"\;",$data[8],100;
		foreach my $info1(@info)
		{
			my($ky,$vl) = split'\=',$info1;
			if($ky eq 'transcript_id')
			{
				$mrnaId = $vl;
				$exonNo{$mrnaId}++;
			}
		}
		if(defined($exonNo{$mrnaId}))
		{
			my $exonNm = $exonNo{$mrnaId};
			if($exonNm>1)
			{
				my $intronStart = $exonEnds{$mrnaId}{($exonNm-1)}+1;
				my $intronEnd = $exonStart-1;
				my $intronLength = 0 ;
				if($data[6] eq '+')
				{
					$intronLength = ($intronEnd-$intronStart)+1;
				}
				elsif($data[6] eq '-')
				{
					$intronLength = ($intronStart-$intronEnd)+1;
				}
				$intronNo{$mrnaId}++;
				my $intronNm = $intronNo{$mrnaId};
				$totalIntrons++;
				$TotalIntronLen+=$intronLength;
				print OUT5 "$mrnaId\t$intronNm\t$intronStart\t$intronEnd\t$intronLength\n";
			}
			$exonEnds{$mrnaId}{$exonNm} = $exonEnd;
			$totalExons++;
			$TotalExonLen+=$exonLength;
			print OUT4 "$mrnaId\t$exonNm\t$exonStart\t$exonEnd\t$exonLength\n";
		}
	}
	elsif($data[2] eq 'CDS')
	{
		my @info = split"\;",$data[8],100;
		my $proteinId = 'NA';
		my $Id = 'Na';
		foreach my $info1(@info)
		{
			my($ky,$vl) = split'\=',$info1;
			if($ky eq 'protein_id')
			{
				$proteinId = $vl;
			}
			if($ky eq 'Parent')
			{
				$Id = $vl;
			}
		}
		my $mrnaId = $mRNAMap{$Id};
		$ProteinInfo{$mrnaId} = $proteinId;
	}
}
close IN1;
close OUT2;
close OUT4;
close OUT5;

# Parsing Protein Table File #
open IN2,"<$infile2" or die "Can't open $infile2 for reading\n";
while(<IN2>)
{
	chomp;
	next if(/^\s*$/);
	next if(/^\s*\#+.*/);
	my @data = split"\t",$_,999;
	my $proteinLength = $data[8];
	my $proteinId = $data[7];
	my $proteinName = $data[9];
	$ProteinDetails{$proteinId} = "$proteinId\t$proteinLength\t$proteinName";
}
close IN2;


open OUT3,">$resulfile3" or die "Can't open $resulfile3 for writing\n";
open OUT6,">$resulfile6" or die "Can't open $resulfile6 for writing\n";
foreach my $mNo(sort{$a<=>$b;} keys %mrnaInfo)
{
	foreach my $mrnaid(keys %{$mrnaInfo{$mNo}})
	{
		$TotalmRNA++;
		my $proteinid = 'NA';
		if(defined($ProteinInfo{$mrnaid}))
		{
			$proteinid = $ProteinInfo{$mrnaid};
			my ($pid,$pLen,$pDesc) = ('NA','NA','NA');
			if(defined($ProteinDetails{$proteinid}))
			{
				($pid,$pLen,$pDesc) = split"\t",$ProteinDetails{$proteinid};
			}
			$TotalProteins++;
			print OUT6 "$mrnaid\t$pid\t$pLen\t$pDesc\n";			
		}
		print OUT3 "".$mrnaInfo{$mNo}{$mrnaid}."\t$proteinid\n";
	}
}
close OUT3;
close OUT6;

open OUT1,">$resulfile1" or die "Can't open $resulfile1 for writing\n";
my $averageExonLength = int($TotalExonLen/$totalExons);
my $averageIntronLength = int($TotalIntronLen/$totalIntrons);
my $averageExonsPerGene = sprintf("%0.2f",($totalExons/$TotalGenes));
my $averageIntronsPerGene = sprintf("%0.2f",($totalIntrons/$TotalGenes));
print OUT1 "$taxid\t$organism\t$commonName\t$TotalGenes\t$TotalmRNA\t$TotalProteins\t$totalExons\t$totalIntrons\t$averageExonLength\t$averageIntronLength\t$averageExonsPerGene\t$averageIntronsPerGene\n";
close OUT1;
