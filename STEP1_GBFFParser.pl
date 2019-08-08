#BEGIN#######################################################
# Version:   1.0
# Filename:  GBKParser.pl
# Author: Bipin Balan
# Date: 14/09/2017
# Purpose: Parsing the Gene information from GBK Files).
# Notes:
#    1.) The Input will be downloaded GBK files.
##############################################################
#       Modification Block
# Number #   Author    Change Description.      Date     Version
#
#########################################################END#
use strict;
use DBI;
use Bio::SeqIO;
use Bio::Perl;
use Bio::SeqFeature::Tools::Unflattener;
use Getopt::Long;

####### PARAMETERS #######################

# Gene2Refseq File path #
my $Gene2RefSeq = "";
my $GBKFile = "";
GetOptions ("gene2refseq=r" => \$Gene2RefSeq,
              "file=f"   => \$GBKFile)
or die("Error in command line arguments\n");
  
my $resultfolder = 'PARSED_INFO';
unless(-d $resultfolder)
{
	mkdir $resultfolder;
}

# Define File Format #
my $format = 'genbank';

# OutputFile #
my $GBKResultFile = "$resultfolder\/walnut_GBK_Result.txt";

###############################################

# Reading gene2refseq file to get NM and corresponding NP ID #
my %NM_TO_NP_Mapping = ();
&TranscriptToProtein($Gene2RefSeq,\%NM_TO_NP_Mapping);

# Opening ResultFile for writing #
open OUT,">$GBKResultFile" or die "Can't open $GBKResultFile for writing\n";
print OUT "ORGANISM\tTAX_ID\tSTRAND\tGeneID\tGeneSymbol\tGeneStart\tGeneEnd\tHNGC_ID\tOMIM_ID\tTRANSCRIPT_ID\tmRNA_CORDINATES\tmRNA_SEQUENCE\tPROTEIN_ID\tCCDS_ID\tCDS_CORDINATES\tCDS_SEQUENCE\t5\'UTR\t3\'UTR\n";


# Creating BioPerl SeqIO object #
my $seqin = Bio::SeqIO->new(-file => $GBKFile,-format => $format);

# Initialization Block #
my $organism = '-';
my $TaxID = '-';
my %GeneInfo = ();
my $Strand = '';
my $GeneID = '-';
my $GeneSymbol = '-';
my $HNGC_ID = '-';
my $OMIM_ID = '-';
my $GeneStart = '-';
my $GeneEnd = '-';
my $mRNACord = '';
my $mRNASeq = '';
my %mRNAInfo = ();
my $NM_ID = '-';
my $CDSCord = '';
my $CDSSeq = '';
my %CDSInfo = ();
my $NP_ID = '-';
my $CCDS_ID = '-';

while (my $seqio = $seqin->next_seq)
{
	foreach my $SeqFeatures($seqio->get_all_SeqFeatures())
	{
		if($SeqFeatures->primary_tag eq 'source')
		{
			foreach my $feature($SeqFeatures->get_all_tags)
			{
				if($feature eq 'organism') # Fetching Organism Name #
				{
					$organism = ($SeqFeatures->get_tag_values($feature))[0];
				}
				elsif($feature eq 'db_xref')
				{
					foreach my $value ($SeqFeatures->get_tag_values($feature))
					{
						if($value=~/taxon:(\d+)\s*/) # Fetching Taxonomy ID #
						{
							$TaxID = $1;
						}
					}
				}
			}
		}
		elsif($SeqFeatures->primary_tag eq 'gene') # Fetching Gene information under gene tag #
		{
			# Printing O/P of the first gene #
			if(scalar(keys %GeneInfo)>0)
			{
				foreach my $GeneIDs(keys %GeneInfo)
				{
					if(scalar(keys %mRNAInfo)>0)
					{
						foreach my $mRNAIDs(keys %{$mRNAInfo{$GeneIDs}})
						{
							my @CDSinfo = ();
							my @geneinfo = $GeneInfo{$GeneIDs};
							my @mRNAinfo = @{$mRNAInfo{$GeneIDs}{$mRNAIDs}};
							my $UTR_5 = '-';
							my $UTR_3 = '-';
							if(defined($CDSInfo{$GeneIDs}{$NM_TO_NP_Mapping{$mRNAIDs}}))
							{
								@CDSinfo = @{$CDSInfo{$GeneIDs}{$NM_TO_NP_Mapping{$mRNAIDs}}};
								my $mrnaStart = (split'\.\.',(split'\,',$mRNAinfo[1])[0])[0];
								my $cdsStart = (split'\.\.',(split'\,',$CDSinfo[2])[0])[0];
								$UTR_5 = ($mrnaStart != $cdsStart)?"$mrnaStart\.\.$cdsStart":'-';
								my $mrnaEnd = (split'\.\.',(split'\,',$mRNAinfo[1])[-1])[1];
								my $cdsEnd = (split'\.\.',(split'\,',$CDSinfo[2])[-1])[1];
								$UTR_3 = ($mrnaEnd != $cdsEnd)?"$cdsEnd\.\.$mrnaEnd":'-';
							}
							else
							{
								@CDSinfo = ('-','-','-','-');
							}
							my $GeneData = join "\t",@geneinfo,@mRNAinfo,@CDSinfo,$UTR_5,$UTR_3;
							print OUT "$GeneData\n";
						}
					}
					elsif(scalar(keys %CDSInfo)>0)
					{
						foreach my $CDSIDs(keys %{$CDSInfo{$GeneIDs}})
						{
							my @geneinfo = $GeneInfo{$GeneIDs};
							my @mRNAinfo = ('-','-','-');
							my @CDSinfo = @{$CDSInfo{$GeneIDs}{$CDSIDs}};
							my $UTR_5 = '-';
							my $UTR_3 = '-';
							my $GeneData = join "\t",@geneinfo,@mRNAinfo,@CDSinfo,$UTR_5,$UTR_3;
							print OUT "$GeneData\n";
						}
					}
				}
			}

			# Re-Initialization of variables #
			%GeneInfo = ();
			$Strand = '';
			$GeneID = '-';
			$GeneSymbol = '-';
			$HNGC_ID = '-';
			$OMIM_ID = '-';
			$GeneStart = '-';
			$GeneEnd = '-';
			$mRNACord = '';
			$mRNASeq = '';
			%mRNAInfo = ();
			$NM_ID = '-';
			$CDSCord = '';
			$CDSSeq = '';
			%CDSInfo = ();
			$NP_ID = '-';
			$CCDS_ID = '-';

			# Getting strand information. (-1 means negative strand and 1 mens positive strand) #
			$Strand = $SeqFeatures->strand();

			# Getting Gene Start and Gene End #
			$GeneStart = $SeqFeatures->start;
			$GeneEnd = $SeqFeatures->end;

			foreach my $feature($SeqFeatures->get_all_tags)
			{
				if($feature eq 'gene')
				{
					$GeneSymbol = ($SeqFeatures->get_tag_values($feature))[0];
				}
				elsif($feature eq 'db_xref')
				{
					foreach my $value ($SeqFeatures->get_tag_values($feature))
					{
						if($value=~/GeneID:(\d+)\s*/)
						{
							$GeneID = $1; # Gene ID #
						}
						elsif($value=~/HGNC:(\d+)\s*/) # HGNC ID #
						{
							$HNGC_ID = $1;
						}
						elsif($value=~/MIM:(\d+)\s*/)
						{
							$OMIM_ID = $1; # OMIM ID #
						}
					}
				}
			}
			$GeneInfo{$GeneID} = "$organism\t$TaxID\t$Strand\t$GeneID\t$GeneSymbol\t$GeneStart\t$GeneEnd\t$HNGC_ID\t$OMIM_ID";
		}
		elsif($SeqFeatures->primary_tag eq 'mRNA') # Fetching mRNA information under mRNA tag #
		{
			$mRNACord = '';
			$mRNASeq = '';

			foreach my $feature($SeqFeatures->get_all_tags)
			{
				if($feature eq 'transcript_id')
				{
					$NM_ID = ($SeqFeatures->get_tag_values($feature))[0];
				}
			}

			if ( $SeqFeatures->location->isa('Bio::Location::SplitLocationI')) # If multiple cordinates #
			{
				for my $location ( $SeqFeatures->location->sub_Location )
				{
					if($mRNACord eq '')
					{
						$mRNACord = $location->start."\.\.".$location->end;
					}
					else
					{
						$mRNACord .= "\,".$location->start."\.\.".$location->end;
					}
					if($Strand == -1)
					{
						my $t = $seqio->subseq($location->start,$location->end);
						$t=~tr/ATGC/TACG/;
						$mRNASeq .= $t;
					}
					else
					{
						$mRNASeq .= $seqio->subseq($location->start,$location->end);
					}
				}
				if($Strand == -1)
				{
					$mRNASeq = reverse($mRNASeq);
				}
			}
			else # Single cordinate #
			{
				$mRNACord.=$SeqFeatures->start."\.\.".$SeqFeatures->end;
				$mRNASeq = $seqio->subseq($SeqFeatures->start,$SeqFeatures->end);
				if($Strand == -1)
				{
					$mRNASeq=~tr/ATGC/TACG/;
					$mRNASeq = reverse($mRNASeq);
				}					
			}
			@{$mRNAInfo{$GeneID}{$NM_ID}} = ($NM_ID,$mRNACord,$mRNASeq);
		}
		elsif($SeqFeatures->primary_tag eq 'CDS') # Fetching CDS information under CDS tag #
		{
			$CDSCord = '';
			$CDSSeq = '';
			$CCDS_ID = '-';

			foreach my $feature($SeqFeatures->get_all_tags)
			{
				if($feature eq 'protein_id')
				{
					$NP_ID = ($SeqFeatures->get_tag_values($feature))[0];
				}
				elsif($feature eq 'db_xref')
				{
					foreach my $value ($SeqFeatures->get_tag_values($feature))
					{
						if($value=~/CCDS:(\S+)\s*/)
						{
							$CCDS_ID = $1; # CCDS ID #
						}
					}
				}
			}
			if ( $SeqFeatures->location->isa('Bio::Location::SplitLocationI')) # If multiple cordinates #
			{
				for my $location ( $SeqFeatures->location->sub_Location )
				{
					if($CDSCord eq '')
					{
						$CDSCord = $location->start."\.\.".$location->end;
					}
					else
					{
						$CDSCord .= "\,".$location->start."\.\.".$location->end;
					}
					$CDSSeq .= $seqio->subseq($location->start,$location->end);
				}
			}
			else # Single cordinate #
			{
				$CDSCord.=$SeqFeatures->start."\.\.".$SeqFeatures->end;
				$CDSSeq = $seqio->subseq($SeqFeatures->start,$SeqFeatures->end);
			}
			@{$CDSInfo{$GeneID}{$NP_ID}} = ($NP_ID,$CCDS_ID,$CDSCord,$CDSSeq);
		}
	}

	# Printing O/P #
	if(scalar(keys %GeneInfo)>0)
	{
		foreach my $GeneIDs(keys %GeneInfo)
		{
			if(scalar(keys %mRNAInfo)>0)
			{
				foreach my $mRNAIDs(keys %{$mRNAInfo{$GeneIDs}})
				{
					my @CDSinfo = ();
					my @geneinfo = $GeneInfo{$GeneIDs};
					my @mRNAinfo = @{$mRNAInfo{$GeneIDs}{$mRNAIDs}};
					my $UTR_5 = '-';
					my $UTR_3 = '-';
					if(defined($CDSInfo{$GeneIDs}{$NM_TO_NP_Mapping{$mRNAIDs}}))
					{
						@CDSinfo = @{$CDSInfo{$GeneIDs}{$NM_TO_NP_Mapping{$mRNAIDs}}};
						my $mrnaStart = (split'\.\.',(split'\,',$mRNAinfo[1])[0])[0];
						my $cdsStart = (split'\.\.',(split'\,',$CDSinfo[2])[0])[0];
						$UTR_5 = ($mrnaStart != $cdsStart)?"$mrnaStart\.\.$cdsStart":'-';
						my $mrnaEnd = (split'\.\.',(split'\,',$mRNAinfo[1])[-1])[1];
						my $cdsEnd = (split'\.\.',(split'\,',$CDSinfo[2])[-1])[1];
						$UTR_3 = ($mrnaEnd != $cdsEnd)?"$cdsEnd\.\.$mrnaEnd":'-';
					}
					else
					{
						@CDSinfo = ('-','-','-','-');
					}
					my $GeneData = join "\t",@geneinfo,@mRNAinfo,@CDSinfo,$UTR_5,$UTR_3;
					print OUT "$GeneData\n";
				}
			}
			elsif(scalar(keys %CDSInfo)>0)
			{
				foreach my $CDSIDs(keys %{$CDSInfo{$GeneIDs}})
				{
					my @geneinfo = $GeneInfo{$GeneIDs};
					my @mRNAinfo = ('-','-','-');
					my @CDSinfo = @{$CDSInfo{$GeneIDs}{$CDSIDs}};
					my $UTR_5 = '-';
					my $UTR_3 = '-';
					my $GeneData = join "\t",@geneinfo,@mRNAinfo,@CDSinfo,$UTR_5,$UTR_3;
					print OUT "$GeneData\n";
				}
			}
		}
	}

	# Re-Initialization of variables #
	%GeneInfo = ();
	$Strand = '';
	$GeneID = '-';
	$GeneSymbol = '-';
	$HNGC_ID = '-';
	$OMIM_ID = '-';
	$GeneStart = '-';
	$GeneEnd = '-';
	$mRNACord = '';
	$mRNASeq = '';
	%mRNAInfo = ();
	$NM_ID = '-';
	$CDSCord = '';
	$CDSSeq = '';
	%CDSInfo = ();
	$NP_ID = '-';
	$CCDS_ID = '-';
}
close OUT; # Closing File handler #
exit;

sub TranscriptToProtein
{
	my $Gene2RefSeqFile = shift;
	my $RefNM_TO_NP_Mapping = shift;
	
	# Opening file and reading #
	open IN,"<$Gene2RefSeqFile" or die "Can't open $Gene2RefSeqFile for reading\n";
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/);
		next if(/^\s*\#Format\:.*/);
		my @data = split"\t",$_;
		if($data[0]==51240)
		{
			$$RefNM_TO_NP_Mapping{$data[3]} = $data[5];
		}
	}
	close IN;
}
