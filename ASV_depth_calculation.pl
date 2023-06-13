#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Take ASV taxonomy file, ASV counts file, and ASV_depthOfTaxonomicInference file to calculate stats on change in power over time.
#Compare to full run Taxonomy file for accuracy.
#
#Stats:
#    %ASVs assigned to each taxonomic level (% of total with depth of inference valid to that level)
#    %Summed study rel abundance assigned to each taxonomic level
#    Accuracy - %ASV assignments agree with taxonomic structure of full run file.

# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("t:c:d:o:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-t = ASV taxonomy file for date test\n";
        print "-c = ASV counts file\n";
        print "-d = ASV_depthOfTaxonomicInference file\n";
        print "-o = Original ASV taxonomy file for accuracy comparison\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %CURRENT;
my %REFERENCE;

my @hierarchy = qw(kingdom phylum class order family genus species);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
open(IN, "<$options{t}") or die "\n\nNADA $options{t} file!!\n\n";
my @data = <IN>; close(IN); shift(@data);
foreach my $line (@data)
    {   chomp($line);
        my @split_line = split('\t', $line);
        my $asv = $split_line[0];
        my $kingdom = $split_line[1];
        my $phylum = $split_line[2];
        my $class = $split_line[3];
        my $order = $split_line[4];
        my $family = $split_line[5];
        my $genus = $split_line[6];
        my $species = $split_line[7];
        if ($asv =~ m/ASV_/)
            {   $CURRENT{$asv}{'kingdom'} = $kingdom;
                $CURRENT{$asv}{'phylum'} = $phylum;
                $CURRENT{$asv}{'class'} = $class;
                $CURRENT{$asv}{'order'} = $order;
                $CURRENT{$asv}{'family'} = $family;
                $CURRENT{$asv}{'genus'} = $genus;
                $CURRENT{$asv}{'species'} = $species;
            }
        my $deepestHierarchy;
        my $deepestSwitch = "TRUE";
        foreach my $i (0..$#hierarchy)
            {   if ($CURRENT{$asv}{$hierarchy[$i]} eq "NA")
                    {   $deepestHierarchy = $hierarchy[$i - 1];
                        $CURRENT{$asv}{'deepest'} = $deepestHierarchy;
                        $CURRENT{$asv}{'deepest_element'} = $i - 1;
                        #DEEPEST = Furthest to the right
                        $deepestSwitch = "FALSE";
                        last;
                    }
            }
        if ($deepestSwitch eq "TRUE")
            {   $CURRENT{$asv}{'deepest'} = "species";
                $CURRENT{$asv}{'deepest_element'} = 6;
            }
    }
    
open(IN2, "<$options{o}") or die "\n\nNADA $options{o} file!!\n\n";
my @data2 = <IN2>; close(IN2); shift(@data2);
foreach my $line (@data2)
    {   chomp($line);
        my @split_line = split('\t', $line);
        my $asv = $split_line[0];
        my $kingdom = $split_line[1];
        my $phylum = $split_line[2];
        my $class = $split_line[3];
        my $order = $split_line[4];
        my $family = $split_line[5];
        my $genus = $split_line[6];
        my $species = $split_line[7];
        if ($asv =~ m/ASV_/)
            {   $REFERENCE{$asv}{'kingdom'} = $kingdom;
                $REFERENCE{$asv}{'phylum'} = $phylum;
                $REFERENCE{$asv}{'class'} = $class;
                $REFERENCE{$asv}{'order'} = $order;
                $REFERENCE{$asv}{'family'} = $family;
                $REFERENCE{$asv}{'genus'} = $genus;
                $REFERENCE{$asv}{'species'} = $species;
            }
        my $deepestHierarchy;
        my $deepestSwitch = "TRUE";
        foreach my $i (0..$#hierarchy)
            {   if ($REFERENCE{$asv}{$hierarchy[$i]} eq "NA")
                    {   $deepestHierarchy = $hierarchy[$i - 1];
                        $REFERENCE{$asv}{'deepest'} = $deepestHierarchy;
                        #DEEPEST = Furthest to the right
                        $deepestSwitch = "FALSE";
                        last;
                    }
            }
        if ($deepestSwitch eq "TRUE")
            {   $REFERENCE{$asv}{'deepest'} = "species";
            }
    }
    
    
open(IN3, "<$options{d}") or die "\n\nNADA $options{d} file!!\n\n";
my @data3 = <IN3>; close(IN3); shift(@data3);
foreach my $line (@data3)
    {   chomp($line);
        my @split_line = split('\t', $line);
        my $asv = $split_line[0];
        my $lowestConfidence = lc($split_line[2]);
        $CURRENT{$asv}{'shallowestConfidenceFromHits'} = $lowestConfidence;
        foreach my $j (0..$#hierarchy)
            {   if ($hierarchy[$j] eq $lowestConfidence)
                    {   $CURRENT{$asv}{'shallowestConfidenceFromHits_element'} = $j;
                        #SHALLOWEST = FURTHEST FROM THE RIGHT
                    }
            }
        
        
        
    }
    
open(IN4, "<$options{c}") or die "\n\nNADA $options{c} file!!\n\n";
my @data4 = <IN4>; close(IN4); shift(@data4);
foreach my $line (@data4)
    {   chomp($line);
        my @split_line = split('\t', $line);
        my $asv = $split_line[0];
        my $totalCount = 0;
        foreach my $i (1..$#split_line)
            {   $totalCount += $split_line[$i];
            }
        $CURRENT{$asv}{'countSum'} = $totalCount;
    }
    
#COMPARE CURRENT AND REFERENCE
open(ACCURACY, ">ASV_accuracy.txt");
foreach my $i (sort keys %CURRENT)
    {   if (exists $CURRENT{$i}{'shallowestConfidenceFromHits_element'})
            {   my $currentLow = $CURRENT{$i}{'deepest'};
                my $currentLowElement = $CURRENT{$i}{'deepest_element'};
                my $confidence = $CURRENT{$i}{'shallowestConfidenceFromHits_element'};
                my $currentCompareTaxa;
                if ($currentLowElement >= $confidence)
                    {   $currentCompareTaxa = $CURRENT{$i}{$currentLow};
                        if ($currentCompareTaxa eq $REFERENCE{$i}{$currentLow})
                            {   $CURRENT{$i}{'accurate'} = "TRUE";
                                print ACCURACY "$i\tTRUE\n";
                            }
                        else
                            {   $CURRENT{$i}{'accurate'} = "FALSE";
                                print ACCURACY "$i\tFALSE\n";
                            }
                    }
                else
                    {   $CURRENT{$i}{'accurate'} = "NA";
                        print ACCURACY "$i\tNA\n";
                    }
            }
        else
            {   $CURRENT{$i}{'accurate'} = "NA";
                print ACCURACY "$i\tNA\n";
            }
    }
close(ACCURACY);
    
    
#ANALYSIS/PRINT TO STDOUT
my @species;
my @genus;
my @family;
my @order;
my @class;

my $totalASVs = 0;
my $AccurateASVs = 0;
foreach my $i (sort keys %CURRENT)
    {   if ($CURRENT{$i}{'accurate'} eq "TRUE" || $CURRENT{$i}{'accurate'} eq "FALSE")
            {   $totalASVs += 1;
            }
        if ($CURRENT{$i}{'accurate'} eq "TRUE")
            {   $AccurateASVs += 1;
            }
        if (exists $CURRENT{$i}{'shallowestConfidenceFromHits'})
            {   my $dep = $CURRENT{$i}{'shallowestConfidenceFromHits'};
                if ($dep eq "species")
                    {   push(@species, $i);
                    }
                if ($dep eq "genus")
                    {   push(@species, $i);
                        push(@genus, $i);
                    }
                if ($dep eq "family")
                    {   push(@species, $i);
                        push(@genus, $i);
                        push(@family, $i);
                    }
                if ($dep eq "order")
                    {   push(@species, $i);
                        push(@genus, $i);
                        push(@family, $i);
                        push(@order, $i);
                    }
                if ($dep eq "class")
                    {   push(@species, $i);
                        push(@genus, $i);
                        push(@family, $i);
                        push(@order, $i);
                        push(@class, $i);
                    }
                if ($dep eq "phylum")
                    {   push(@species, $i);
                        push(@genus, $i);
                        push(@family, $i);
                        push(@order, $i);
                        push(@class, $i);
                    }
                if ($dep eq "kingdom")
                    {   push(@species, $i);
                        push(@genus, $i);
                        push(@family, $i);
                        push(@order, $i);
                        push(@class, $i);
                    }
            }
    }

my $valid_total_CouldBespecies = scalar(@species);
my $valid_total_CouldBegenus = scalar(@genus);
my $valid_total_CouldBefamily = scalar(@family);
my $valid_total_CouldBeorder = scalar(@order);
my $valid_total_CouldBeclass = scalar(@class);

my $sum_counts_CouldBespecies = 0;
my $sum_counts_CouldBegenus = 0;
my $sum_counts_CouldBefamily = 0;
my $sum_counts_CouldBeorder = 0;
my $sum_counts_CouldBeclass = 0;
foreach my $i (@species)
    {   $sum_counts_CouldBespecies += $CURRENT{$i}{'countSum'};   
    }
foreach my $i (@genus)
    {   $sum_counts_CouldBegenus += $CURRENT{$i}{'countSum'};   
    }
foreach my $i (@family)
    {   $sum_counts_CouldBefamily += $CURRENT{$i}{'countSum'};   
    }
foreach my $i (@order)
    {   $sum_counts_CouldBeorder += $CURRENT{$i}{'countSum'};   
    }
foreach my $i (@class)
    {   $sum_counts_CouldBeclass += $CURRENT{$i}{'countSum'};   
    }


my $valid_asvs_actually_species = 0;
my $valid_asvs_actually_genus = 0;
my $valid_asvs_actually_family = 0;
my $valid_asvs_actually_order = 0;
my $valid_asvs_actually_class = 0;
my $sum_counts_Actuallyspecies = 0;
my $sum_counts_Actuallygenus = 0;
my $sum_counts_Actuallyfamily = 0;
my $sum_counts_Actuallyorder = 0;
my $sum_counts_Actuallyclass = 0;

foreach my $i (sort keys %CURRENT)
    {   my $validspecies = "FALSE";
        my $validgenus = "FALSE";
        my $validfamily = "FALSE";
        my $validorder = "FALSE";
        my $validclass = "FALSE";
        foreach my $j (@species)
            {   if ($i eq $j)
                    {   $validspecies = "TRUE";
                    }
            }
        foreach my $j (@genus)
            {   if ($i eq $j)
                    {   $validgenus = "TRUE";
                    }
            }
        foreach my $j (@family)
            {   if ($i eq $j)
                    {   $validfamily = "TRUE";
                    }
            }
        foreach my $j (@order)
            {   if ($i eq $j)
                    {   $validorder = "TRUE";
                    }
            }
        foreach my $j (@class)
            {   if ($i eq $j)
                    {   $validclass = "TRUE";
                    }
            }
        if ($CURRENT{$i}{'species'} ne "NA" && $validspecies eq "TRUE")
            {   $valid_asvs_actually_species += 1;
                $sum_counts_Actuallyspecies += $CURRENT{$i}{'countSum'};
            }
        if ($CURRENT{$i}{'genus'} ne "NA" && $validgenus eq "TRUE")
            {   $valid_asvs_actually_genus += 1;
                $sum_counts_Actuallygenus += $CURRENT{$i}{'countSum'};
            }
        if ($CURRENT{$i}{'family'} ne "NA" && $validfamily eq "TRUE")
            {   $valid_asvs_actually_family += 1;
                $sum_counts_Actuallyfamily += $CURRENT{$i}{'countSum'};
            }
        if ($CURRENT{$i}{'order'} ne "NA" && $validorder eq "TRUE")
            {   $valid_asvs_actually_order += 1;
                $sum_counts_Actuallyorder += $CURRENT{$i}{'countSum'};
            }
        if ($CURRENT{$i}{'class'} ne "NA" && $validclass eq "TRUE")
            {   $valid_asvs_actually_class += 1;
                $sum_counts_Actuallyclass += $CURRENT{$i}{'countSum'};
            }
    }
    
print "TotalASVs_species\tTotalASVs_genus\tTotalASVs_family\tTotalASVs_order\tTotalASVs_class\t";
print "%ASVs_species\t%ASVs_genus\t%ASVs_family\t%ASVs_order\t%ASVs_class\t";
print "%Community_species\t%Community_genus\t%Community_family\t%Community_order\t%Community_class\t";
print "%AccurateASVs\tTotalASVsinAccuracyEstimate\n";

print "$valid_total_CouldBespecies\t$valid_total_CouldBegenus\t$valid_total_CouldBefamily\t$valid_total_CouldBeorder\t$valid_total_CouldBeclass\t";

my $round_ASVs_species;
my $round_ASVs_genus;
my $round_ASVs_family;
my $round_ASVs_order;
my $round_ASVs_class;

if ($valid_total_CouldBespecies == 0)
    {$round_ASVs_species = "NA";}
else
    {   my $percASVs_species = (100 * $valid_asvs_actually_species) / $valid_total_CouldBespecies;
        $round_ASVs_species = &ROUND($percASVs_species,2);
    }
if ($valid_total_CouldBegenus == 0)
    {$round_ASVs_genus = "NA";}
else
    {   my $percASVs_genus = (100 * $valid_asvs_actually_genus) / $valid_total_CouldBegenus;
        $round_ASVs_genus = &ROUND($percASVs_genus,2); 
    }
if ($valid_total_CouldBefamily == 0)
    {$round_ASVs_family = "NA";}
else
    {   my $percASVs_family = (100 * $valid_asvs_actually_family) / $valid_total_CouldBefamily;
        $round_ASVs_family = &ROUND($percASVs_family,2);
    }
if ($valid_total_CouldBeorder == 0)
    {$round_ASVs_order = "NA";}
else
    {   my $percASVs_order = (100 * $valid_asvs_actually_order) / $valid_total_CouldBeorder;
        $round_ASVs_order = &ROUND($percASVs_order,2);
    }
if ($valid_total_CouldBeclass == 0)
    {$round_ASVs_class = "NA";}
else
    {   my $percASVs_class = (100 * $valid_asvs_actually_class) / $valid_total_CouldBeclass;
        $round_ASVs_class = &ROUND($percASVs_class,2);
    }

print "$round_ASVs_species\t$round_ASVs_genus\t$round_ASVs_family\t$round_ASVs_order\t$round_ASVs_class\t";


my $round_Comm_species;
my $round_Comm_genus;
my $round_Comm_family;
my $round_Comm_order;
my $round_Comm_class;

if ($sum_counts_CouldBespecies == 0)
    {$round_Comm_species = "NA";}
else
    {   my $percComm_species = (100 * $sum_counts_Actuallyspecies) / $sum_counts_CouldBespecies;
        $round_Comm_species = &ROUND($percComm_species,2);
    }
if ($sum_counts_CouldBegenus == 0)
    {$round_Comm_genus = "NA";}
else
    {   my $percComm_genus = (100 * $sum_counts_Actuallygenus) / $sum_counts_CouldBegenus;
        $round_Comm_genus = &ROUND($percComm_genus,2);
    }
if ($sum_counts_CouldBefamily == 0)
    {$round_Comm_family = "NA";}
else
    {   my $percComm_family = (100 * $sum_counts_Actuallyfamily) / $sum_counts_CouldBefamily;
        $round_Comm_family = &ROUND($percComm_family,2);
    }
if ($sum_counts_CouldBeorder == 0)
    {$round_Comm_order = "NA";}
else
    {   my $percComm_order = (100 * $sum_counts_Actuallyorder) / $sum_counts_CouldBeorder;
        $round_Comm_order = &ROUND($percComm_order,2);
    }
if ($sum_counts_CouldBeclass == 0)
    {$round_Comm_class = "NA";}
else
    {   my $percComm_class = (100 * $sum_counts_Actuallyclass) / $sum_counts_CouldBeclass;
        $round_Comm_class = &ROUND($percComm_class,2);
    }

print "$round_Comm_species\t$round_Comm_genus\t$round_Comm_family\t$round_Comm_order\t$round_Comm_class\t";

my $percAccurate = (100 * $AccurateASVs) / $totalASVs;
my $round_percAccurate = &ROUND($percAccurate,2);

print "$round_percAccurate\t$totalASVs\n";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub ROUND
{	my $num = shift(@_);
	my $decimals = shift(@_);
	my $roundnum = int(($num * 10**$decimals) + 0.5)/(10**$decimals);
	return $roundnum;
}

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
