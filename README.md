# Filtering BLAST data by release date

Written by Sean McAllister, v.1 01/2023

## Purpose

The purpose of this repository is to encapsulate the methods and scripts used to reformat tab-delimited (-outfmt '6 qseqid pident length staxids sacc') BLAST outfiles to exclude matches released by NCBI after a set date, and to perform this iteratively to assess the ability of the NCBI nt reference database to resolve metabarcoding Amplicon Sequence Variants (ASVs) to Class, Order, Family, Genus, and Species taxonomic levels. Uses REVAMP to assign taxonomy (available here: URL).

## Files

1. Example workflow. 
2. btab_date_manipulations.pl – Remove blast hits in btab file based on reference database.
3. reformat_blast.R – Reformat blast output to choose best hit longer than user-specified cutoff.
4. ASV_depth_calculation.pl – Calculates depths of taxonomic inference.
5. reference_databases – Example pre-calculated reference databases from prior OME runs.

## Workflow
The workflow for this processing progresses through the following steps:
1. Download date information from NCBI for all blast hit accessions.
2. Run btab_date_manipulations.pl iteratively to first develop the Reference Date/Accession Database for gene of interest. When complete, no MissingDate_Accessions.txt will be created for follow up.
3. Reformat the btab file to choose the best blast hits.
4. Run REVAMP, copying off important files.
5. Run ASV_depth_calculation.pl script to calculate Total ASVs to taxa, Percent ASVs to taxa, Percent of overall community to taxa, and Percent Accuracy to current taxonomic assignments.
6. Clean up final files.

## btab_date_manipulations.pl
Takes blastn ```-outfmt 6``` btab data file and produces new file with all records submitted to GenBank after a particular date excluded. Deciphers date information from GenBank records. Prints to file accessions missing date information for follow up. Prints to file all accession/date decisions for a reference database (so future runs do not need additional compute time for the same accession).

Assesses from the btab file the maximum depth appropriate for taxonomic inference for each ASV. This is assigned based on user-supplied taxonomic cutoffs. The maximum depth for inference is defined as the taxonomic level one up from the level inclusive of the lowest percent identity (PID) BLAST assignment. Thus, of a list of (e.g. 4000) BLAST hits with the lowest PID hit at the ORDER level, the depth of taxonomic inference is set to FAMILY, assuming all possible FAMILY level hits are captured in the BLAST data. This is important when BLAST hits are deleted from the btab file, ensuring that all results from a particular taxonomic depth are included before assessment over time. Otherwise, a lack of taxonomic assignment may be the result of a high number of deleted over-represented sequences with no next-best hits to replace them because they were not recorded in the original btab file. 

Use the same taxonomy cutoff string ```-t``` as used in REVAMP. Example for protein encoding genes: "95,92,87,77,67,60". Example for ribosomal RNA genes: "97,95,90,80,70,60". Cutoffs define the boundaries of Species, Genus, Family, Order, Class, Phylum, and unknown assignments.

### Arguments
```
-i = Input btab BLASTn file (-outfmt '6 qseqid pident length staxids sacc')
-m = Month cuttoff setting (i.e. 1 or 6 or 12)
-y = Year cuttoff setting (i.e. 2014)
-t = Taxonomy cutoffs string (i.e. 95,92,87,77,67,60)
-d = Accession and date file (optional if -r provided)
-o = Out directory
-r = Reference file with Acc\tYear\tNumericMonth data (optional)
-h = This help message
```

```-d``` file is produced through ```grep "Submitted (\|ACCESSION\|JOURNAL\|LOCUS"``` from each accession numbers GenBank record.

## reformat_blast.R
Takes BLASTn btab file and simplifies output selecting the best PID greater than a user defined length (best if ~90% of the marker gene length).

### Arguments (positional)
```
POSITION 1 = path to btab file
POSITION 2 = path to out directory
POSITION 3 = numeric denoting the bp length cutoff for consideration
```

## ASV_depth_calculation.pl
Takes REVAMP ASV taxonomy file, ASV counts file, and ASV depth of taxonomic inference file and calculates statistics on Total ASVs assigned to each taxa (Species to Class), Percent of ASVs assigned to each taxa (Species to Class), Percent of overal community (based on ASV relative abundance) to each taxa (Species to Class), and Percent accuracy to current taxonomic assignments (with count of total ASVs used for assessment).

### Arguments
```
-t = ASV taxonomy file for date test (from date run of REVAMP)
-c = ASV counts file (from original run of REVAMP)
-d = ASV_depthOfTaxonomicInference file (from btab_date_manipulations.pl run)
-o = Original ASV taxonomy file for accuracy comparison (from original run of REVAMP)
-h = This help message
```

## Reference Databases
Accession and date reference databases for marker genes already run through this pipeline. Currently: ```Kelly 16S```, ```Machida 18S```, ```Leray COI```, ```Miya MiFish```.

## Dependencies

See REVAMP (URL) documentation for dependencies from that pipeline.

These scripts/workflow require:
1. REVAMP (URL)
2. ```efetch``` (https://www.ncbi.nlm.nih.gov/books/NBK179288/)
3. ```perl``` ```List::MoreUtils```
4. ```Rscript```
5. ```R``` package ```dplyr```

## Outputs
The primary output of this pipeline is a tab-delimited file that requires some cleanup. Each YYYY_MM combiation prints out the following for assessment of the depth of taxonomic assignment and accuracy (compared to current taxonomic assignments) over time. A column might have an "NA" if there are no ASVs where depth can be inferred to that level.

1. ```TotalASVs_species``` = Total number of ASVs that can be assigned to Species.
2. ```TotalASVs_genus``` = Total number of ASVs that can be assigned to Genus.
3. ```TotalASVs_family``` = Total number of ASVs that can be assigned to Family.
4. ```TotalASVs_order``` = Total number of ASVs that can be assigned to Order.
5. ```TotalASVs_class``` = Total number of ASVs that can be assigned to Class.
6. ```%ASVs_species``` = Percent of ASVs that can be assigned to Species.
7. ```%ASVs_genus``` = Percent of ASVs that can be assigned to Genus.
8. ```%ASVs_family``` = Percent of ASVs that can be assigned to Family.
9. ```%ASVs_order``` =  Percent of ASVs that can be assigned to Order.
10. ```%ASVs_class``` = Percent of ASVs that can be assigned to Class.
11. ```%Community_species``` = Percent of the overall community that can be assigned to Species.
12. ```%Community_genus``` = Percent of the overall community that can be assigned to Genus.
13. ```%Community_family``` = Percent of the overall community that can be assigned to Family.
14. ```%Community_order``` = Percent of the overall community that can be assigned to Order.
15. ```%Community_class``` = Percent of the overall community that can be assigned to Class.
16. ```%AccurateASVs``` = Percent of total ASVs that match with current taxonomic assignments.
17. ```TotalASVsinAccuracyEstimate``` = Total number of ASVs that are included in the accuracy estimate.

### Legal Disclaimer

*This repository is a software product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce (DOC). All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the DOC or DOC bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the DOC. The DOC seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by the DOC or the United States Government.*
