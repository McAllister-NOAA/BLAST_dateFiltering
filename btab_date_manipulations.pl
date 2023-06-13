#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Take btab blast file from REVAMP (should have a large number of target subjects)
#Assumption: not going so far back that don't have ALL representative sequences from same order (>75% identity)
#May have to run REVAMP with modified taxa confidence thresholds to assign to "Unknown" at appropriate level (or run separately)
#I.E. So you don't run out of sequences in the btab file and assume that means unknown when it is actually just all of the target sequences we printed to file.

#Due to the time it takes to fetch records from GenBank, I've set up a master file in the format Accession\tYear\tNumericMonth.
#This should be loaded first to find new accessions for a run, to limit the time requirements.

#RUN EFETCH BEFORE THIS PROGRAM (10000 accessions max per query seems to be most efficient, though many fail at 3000)
# cat ASV_blastn_nt.btab | cut -f5 | sort | uniq > uniq_accessions.txt
# split -l 10000 uniq_accessions.txt
# for f in x*; do while read line; do printf "%s" "$line," >>
#     accessionList_${f}.txt; done < $f; efetch -db nuccore -id `cat accessionList_${f}.txt` -format gb -chr_start 1 -chr_stop 5 > genbank_record_${f}.txt; done
# for f in genbank_record*; do grep "Submitted (\|ACCESSION\|JOURNAL\|LOCUS" $f >> accession_date_data.txt; done

#RUN AFTER THIS PROGRAM:
#Rscript --vanilla ~/software/scripts/btab_date_manipulations/reformat_blast.R ASV_blastn_nt_2000_6.btab ./ 300 #Positional arguments are 1) btab file, 2) outdir, 3) blast length cutoff
#Move the "ASV_blastn_nt_formatted.txt" file into the /blast_results directory.
#Delete REVAMP progress until after the reformatting, then rerun. Copy off files you want, and run again.

# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("i:m:t:y:d:o:r:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = Input btab BLASTn file\n";
        print "-m = Month cuttoff setting (i.e. 1 or 6 or 12)\n";
        print "-y = Year cuttoff setting (i.e. 2014)\n";
        print "-t = Taxonomy cutoffs string (i.e. 95,92,87,77,67,60)\n";
        print "-d = Accession and date file (optional if -r provided)\n";
        print "-o = Out directory\n";
        print "-r = Reference file with Acc, Year, Numeric Month data (optional)\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %ACCESSION;
my %GENRECORDS;
my %ASV;
my $yearCutoff = $options{y};
my $monthCutoff = $options{m};
my $missingSwitch = "FALSE";

my @months = qw(JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
system("mkdir -p ".$options{o});
#1) Import reference accession to date master database (Accession\tYear\tNumericMonth)
if ($options{r})
    {   open(REF, "<$options{r}") or die "\n\nNADA $options{r} file!!\n\n";
        my @refdata = <REF>; close(REF);
        foreach my $line (@refdata)
            {   chomp($line);
                my @splitLine = split('\t', $line);
                $ACCESSION{$splitLine[0]}{'year'} = $splitLine[1];
                $ACCESSION{$splitLine[0]}{'month'} = $splitLine[2];
            }
    }
my $newDatabase = "NO";

#2) Import accession/date file to %GENRECORDS.
if ($options{d})
    {
open(IN, "<$options{d}") or die "\n\nNADA $options{d} file!!\n\n";
my @data2 = <IN>; close(IN);
my $accession = "";
my $locusAccession = "";
foreach my $line (@data2)
    {   if ($line =~ m/^LOCUS/)
            {   chomp($line);
                $line =~ m/^LOCUS[ ]+([a-zA-Z0-9_]+)[ ]+.+[0-9]+\-([A-Za-z]+)\-([0-9]+)$/;
                $locusAccession = $1;
                my $month = uc($2);
                my $year = $3;
                my $numericmonth;
                my $monthhit = "FALSE";
                foreach my $i (0..11)
                    {   if ($months[$i] eq $month)
                            {   $numericmonth = $i + 1;
                                $monthhit = "TRUE";
                            }
                    }
                if ($monthhit eq "FALSE")
                    {   print "MONTH Abbreviation not found: $locusAccession. Set to JAN.\n";
                        $numericmonth = 1;
                    }
                $GENRECORDS{$locusAccession}{'locus_acc'} = $locusAccession;
                $GENRECORDS{$locusAccession}{'locus_year'} = $year;
                $GENRECORDS{$locusAccession}{'locus_month'} = $numericmonth;
                $GENRECORDS{$locusAccession}{'match_locus_record_acc'} = "FALSE";
            }
        if ($line =~ m/^ACCESSION/)
            {   chomp($line);
                $line =~ m/^ACCESSION[ ]+([a-zA-Z0-9_]+)/;
                $accession = $1;
                $GENRECORDS{$accession}{'acc'} = $accession;
                if (exists $GENRECORDS{$accession}{'locus_acc'} && $accession eq $locusAccession)
                    {   $GENRECORDS{$accession}{'match_locus_record_acc'} = "TRUE";
                    }
                else
                    {   delete $GENRECORDS{$locusAccession};
                    }
            }
        if ($line =~ m/JOURNAL/ || $line =~ m/Submitted \([0-9]+\-[A-Za-z]+\-[0-9]+\)/)
            {   chomp($line);
                my $month;
                my $year;
                my $numericmonth;
                if ($line =~ m/Submitted \([0-9]+\-[A-Za-z]+\-[0-9]+\)/)
                    {   $line =~ m/Submitted \([0-9]+\-([A-Za-z]+)\-([0-9]+)\)/;
                        $month = uc($1);
                        $year = $2;
                        my $monthhit = "FALSE";
                        foreach my $i (0..11)
                            {   if ($months[$i] eq $month)
                                    {   $numericmonth = $i + 1;
                                        $monthhit = "TRUE";
                                    }
                            }
                        if ($monthhit eq "FALSE")
                            {   print "MONTH Abbreviation not found: $accession. Set to JAN.\n";
                                $numericmonth = 1;
                            }
                        if (defined $accession)
                            {   $GENRECORDS{$accession}{'submitted_year'} .= $year.";";
                                $GENRECORDS{$accession}{'submitted_month'} .= $numericmonth.";";
                                $GENRECORDS{$accession}{'preferred'} = "TRUE";
                            }
                    }
                elsif ($line =~ m/\([12][09][0-9][0-9]\)$/)
                    {   $line =~ m/\(([12][09][0-9][0-9])\)$/;
                        $numericmonth = 1;
                        $year = $1;
                        if (defined $accession)
                            {   $GENRECORDS{$accession}{'pub_year'} .= $year.";";
                                $GENRECORDS{$accession}{'pub_month'} .= $numericmonth.";";
                                $GENRECORDS{$accession}{'pubOnly'} = "TRUE";
                            }
                    }
                elsif ($line =~ m/\([12][09][0-9][0-9]\) In press$/)
                    {   $line =~ m/\(([12][09][0-9][0-9])\) In press$/;
                        $numericmonth = 1;
                        $year = $1;
                        if (defined $accession)
                            {   $GENRECORDS{$accession}{'pub_year'} .= $year.";";
                                $GENRECORDS{$accession}{'pub_month'} .= $numericmonth.";";
                                $GENRECORDS{$accession}{'pubOnly'} = "TRUE";
                            }
                    }
                elsif ($line =~ m/Thesis \([12][09][0-9][0-9]\) /)
                    {   $line =~ m/Thesis \(([12][09][0-9][0-9])\) /;
                        $numericmonth = 1;
                        $year = $1;
                        if (defined $accession)
                            {   $GENRECORDS{$accession}{'pub_year'} .= $year.";";
                                $GENRECORDS{$accession}{'pub_month'} .= $numericmonth.";";
                                $GENRECORDS{$accession}{'pubOnly'} = "TRUE";
                            }
                    }
            }
    }
   
#2.5) Make accession/date decisions for %ACCESSION
foreach my $i (sort keys %GENRECORDS)
    {   my $accession = $i;
        unless (exists $ACCESSION{$accession})
            {   if (exists $GENRECORDS{$i}{'preferred'} && $GENRECORDS{$i}{'preferred'} eq "TRUE")
                    {   my $temp_year = $GENRECORDS{$i}{'submitted_year'};
                        my $temp_month = $GENRECORDS{$i}{'submitted_month'};
                        my @arr_year = split(';', $temp_year);
                        my @arr_month = split(';', $temp_month);
                        my $choice_year = 5000;
                        my $choice_element;
                        foreach my $j (0..$#arr_year)
                            {   if ($arr_year[$j] < $choice_year)
                                    {   $choice_year = $arr_year[$j];
                                        $choice_element = $j;
                                    }
                            }
                        my $finalYear = $choice_year;
                        my $finalMonth = $arr_month[$choice_element];
                        $ACCESSION{$i}{'year'} = $finalYear;
                        $ACCESSION{$i}{'month'} = $finalMonth;
                    }
                elsif (exists $GENRECORDS{$i}{'pubOnly'} && $GENRECORDS{$i}{'pubOnly'} eq "TRUE")
                    {   my $temp_year = $GENRECORDS{$i}{'pub_year'};
                        my $temp_month = $GENRECORDS{$i}{'pub_month'};
                        my @arr_year = split(';', $temp_year);
                        my @arr_month = split(';', $temp_month);
                        my $choice_year = 5000;
                        my $choice_element;
                        foreach my $j (0..$#arr_year)
                            {   if ($arr_year[$j] < $choice_year)
                                    {   $choice_year = $arr_year[$j];
                                        $choice_element = $j;
                                    }
                            }
                        my $finalYear = $choice_year;
                        my $finalMonth = $arr_month[$choice_element];
                        $ACCESSION{$i}{'year'} = $finalYear;
                        $ACCESSION{$i}{'month'} = $finalMonth;
                    }
                elsif (exists $GENRECORDS{$i}{'match_locus_record_acc'} && $GENRECORDS{$i}{'match_locus_record_acc'} eq "TRUE")
                    {   $ACCESSION{$i}{'year'} = $GENRECORDS{$i}{'locus_year'};
                        $ACCESSION{$i}{'month'} = $GENRECORDS{$i}{'locus_month'};                       
                    }
                $newDatabase = "YES";
            }
    }
    }

#3) Read and print ASV_blastn_nt.btab with entries deleted that are newer than the designated month/year cutoffs.
#IF ANY ACCESSIONS ARE MISSING DATES: Print out list one per line (so can be rerun or check manually)
open(IN2, "<$options{i}") or die "\n\nNADA $options{i} file!!\n\n";
my @data = <IN2>; close(IN2);
open(MISSING, ">".$options{o}."/MissingDate_Accessions.txt");
open(FINAL, ">".$options{o}."/ASV_blastn_nt_".$yearCutoff."_".$monthCutoff.".btab");
foreach my $line (@data)
    {   chomp($line);
        my @split_line = split('\t', $line);
        my $asv = $split_line[0];
        my $pid = $split_line[1];
        my $acc = $split_line[4];
        #STORE LOWEST PID PER ACCESSION for later lookup
        if (exists $ASV{$asv})
            {   my $compare = $ASV{$asv}{'pid'};
                if ($compare > $pid)
                    {   $ASV{$asv}{'pid'} = $pid;   
                    }
            }
        else
            {   $ASV{$asv}{'pid'} = $pid;
            }
        #PRINT OUT BTAB FILE IF LESS THAN or EQUAL TO CUTOFFS
        if (exists $ACCESSION{$acc})
            {   if ($ACCESSION{$acc}{'year'} < $yearCutoff)
                    {   if ($missingSwitch eq "FALSE")
                            {   print FINAL "$line\n";
                            }
                    }
                if ($ACCESSION{$acc}{'year'} == $yearCutoff)
                    {   if ($ACCESSION{$acc}{'month'} <= $monthCutoff)
                            {   if ($missingSwitch eq "FALSE")
                                    {   print FINAL "$line\n";
                                    }
                            }
                    }
            }
        else
            {   print MISSING "$acc\n";
                if ($missingSwitch eq "FALSE")
                    {   print FINAL "SOME GENBANK RECORDS ARE MISSING\nRUN ABORTED\n";
                    }
                $missingSwitch = "TRUE";
            }
    }
close(MISSING);
close(FINAL);

#4) Assess btab results for depth of taxonomic inference.
open(DEPTH, ">".$options{o}."/ASV_depthOfTaxonomicInference.txt");
print DEPTH "ASV\tLowestPID\tDepthTaxonomicInference\n";
foreach my $i (sort keys %ASV)
    {   my $filterinfo = $options{t};
		my @splitfilter = split(',', $filterinfo);
		if ($ASV{$i}{'pid'} >= $splitfilter[0])
			{	print DEPTH "$i\t$ASV{$i}{'pid'}\tnotDeepEnough\n";	
			}
		if ($ASV{$i}{'pid'} < $splitfilter[0] && $ASV{$i}{'pid'} >= $splitfilter[1])
			{	print DEPTH "$i\t$ASV{$i}{'pid'}\tSPECIES\n";
			}
		if ($ASV{$i}{'pid'} < $splitfilter[1] && $ASV{$i}{'pid'} >= $splitfilter[2])
			{	print DEPTH "$i\t$ASV{$i}{'pid'}\tGENUS\n";
			}
        if ($ASV{$i}{'pid'} < $splitfilter[2] && $ASV{$i}{'pid'} >= $splitfilter[3])
			{	print DEPTH "$i\t$ASV{$i}{'pid'}\tFAMILY\n";
			}
        if ($ASV{$i}{'pid'} < $splitfilter[3] && $ASV{$i}{'pid'} >= $splitfilter[4])
			{	print DEPTH "$i\t$ASV{$i}{'pid'}\tORDER\n";
			}
        if ($ASV{$i}{'pid'} < $splitfilter[4] && $ASV{$i}{'pid'} >= $splitfilter[5])
			{	print DEPTH "$i\t$ASV{$i}{'pid'}\tCLASS\n";
			}
		if ($ASV{$i}{'pid'} < $splitfilter[5])
			{	print DEPTH "$i\t$ASV{$i}{'pid'}\tPHYLUM\n";
			}
    }
close(DEPTH);

#5) Print new master database file
if ($newDatabase eq "YES")
    {   open(REFOUT, ">".$options{o}."/Reference_DateAccession_Database.txt");
        foreach my $i (sort keys %ACCESSION)
            {   print REFOUT "$i\t$ACCESSION{$i}{'year'}\t$ACCESSION{$i}{'month'}\n";
            }
        close(REFOUT);
    }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
