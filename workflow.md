# Workflow for date filtering of BLAST database and run through REVAMP

## Summary
The workflow for this processing progresses through the following steps:
1. Download date information from NCBI for all blast hit accessions.
2. Run btab_date_manipulations.pl iteratively to first develop the Reference Date/Accession Database for gene of interest. When complete, no MissingDate_Accessions.txt will be created for follow up.
3. Reformat the btab file to choose the best blast hits.
4. Run REVAMP, copying off important files.
5. Run ASV_depth_calculation.pl script to calculate Total No. ASVs to taxa, Percent of ASVs to taxa, Percent of Overall Community to taxa, and Percent Accuracy to current taxonomic assignments.
6. Clean up final files.

## Workflow

### Download date information from NCBI for all blast hit accessions.

In the example below, the input file is created by blastn using ```-outfmt '6 qseqid pident length staxids sacc'```.


```
cat ASV_blastn_nt.btab | cut -f5 | sort | uniq > uniq_accessions.txt
split -l 10000 uniq_accessions.txt

for f in x*
  do while read line
    do printf "%s" "$line," >> accessionList_${f}.txt
    done < $f 
  efetch -db nuccore -id `cat accessionList_${f}.txt` -format gb -chr_start 1 -chr_stop 5 > genbank_record_${f}.txt
  done

for f in genbank_record*
  do grep "Submitted (\|ACCESSION\|JOURNAL\|LOCUS" $f >> accession_date_data.txt
  done
```

### Run btab_date_manipulations.pl

In the example below, ```-t``` denotes user-selected cutoffs delineating taxonomic hierarchy. The example shown is one we frequently use for ribosomal RNA gene markers. See README for information on other options.

When run the first time, the ```accession_date_data.txt``` file from the previous step is provided to ```-d```.
```
for f in {1995..2022}
  do for g in {6,12}
    do perl btab_date_manipulations.pl -i ASV_blastn_nt.btab -m $g -y $f -t 97,95,90,80,70,60 -o ${f}_${g}_out -d accession_date_data.txt
    done
  done
```

After the first run, if the file ```MissingDate_Accessions.txt``` is created, then there are some accessions in your btab file that don't exist in the ```accession_date_data.txt``` file. Explore this file to determine why dates were not able to be automatically pulled from the GenBank records. Add these records to the ```Reference_DateAccession_Database.txt``` file in the format: ```Accession\tYear(four digits)\tMonth (numeric; no leading zeros)```. Then run again including the ```-r``` option.
```
for f in {1995..2022}
  do for g in {6,12}
    do perl btab_date_manipulations.pl -i ASV_blastn_nt.btab -m $g -y $f -t 97,95,90,80,70,60 -o ${f}_${g}_out -d accession_date_data.txt -r Reference_DateAccession_Database.txt
    done
  done
```

Once the reference database is complete and capturing all accessions from the btab file, you can provide that database only as input.
```
for f in {1995..2022}
  do for g in {6,12}
    do perl btab_date_manipulations.pl -i ASV_blastn_nt.btab -m $g -y $f -t 97,95,90,80,70,60 -o ${f}_${g}_out -r Reference_DateAccession_Database.txt
    done
  done
```

Once happy with output, clean up duplicate files and remove "MissingDate_Accessions.txt" from first runs.

```
cp 2022_12_out/ASV_depthOfTaxonomicInference.txt ./

for f in {1995..2022}
  do for g in {6,12}
    do rm ${f}_${g}_out/ASV_depthOfTaxonomicInference.txt
    done
  done

for f in {1995..2022}
  do for g in {6,12}
    do rm ${f}_${g}_out/MissingDate_Accessions.txt
    done
  done
```

### Reformat the btab file to choose the best blast hits

Modify with your own paths.

```
for f in {1995..2022}
  do for g in {6,12}
    do Rscript --vanilla reformat_blast.R /use/absolute/PATH/${f}_${g}_out/ASV_blastn_nt_${f}_${g}.btab /use/absolute/PATH/${f}_${g}_out 300
    done
  done
```

### Run REVAMP, copying off important files
Note: It is a good idea to copy off your original REVAMP run before doing these manipulations on the same folder.

```
for f in {1995..2022}
  do for g in {6,12}
    do workingDir=~/Desktop/REVAMP_folder
    
    #Moves formatted blast file from the date formatted folder to the primary REVAMP blast_results directory
    cp ${workingDir}/blast_results/${f}_${g}_out/ASV_blastn_nt_formatted.txt ${workingDir}/blast_results/
    
    #Reforms the progress.txt file to reset for the next date run
    echo "cutadaptFinished=TRUE" > ${workingDir}/progress.txt
    echo "dada2_Finished=TRUE" >> ${workingDir}/progress.txt
    echo "blastFinished=TRUE" >> ${workingDir}/progress.txt
    echo "blastformattingFinished=TRUE" >> ${workingDir}/progress.txt
    
    #Run REVAMP
    revamp.sh -p config_file.txt -f figure_config_file.txt -s samplemetadata.txt -r reads -o OUT -y
    mkdir -p ${workingDir}/blast_results/date_testing/${f}_${g}_out
    cp -r ${workingDir}/processed_tab* ${workingDir}/blast_results/date_testing/${f}_${g}_out/
    cp -r ${workingDir}/ASV2Taxono* ${workingDir}/blast_results/date_testing/${f}_${g}_out/
    cp -r ${workingDir}/Figur* ${workingDir}/blast_results/date_testing/${f}_${g}_out/ #Optional
    done
  done
```

### Run ASV_depth_calculation.pl script to calculate stats
Note that files are coming from three sources:
1. Original REVAMP run for ```-c``` and ```-o```.
2. ```ASV_depthOfTaxonomicInference.txt``` created after completion of ```btab_date_manipulations.pl``` iteration for ```-d```.
3. Specific ASV taxonomy assignment table from each date run for ```-t```.

```
for f in {1995..2022}
  do for g in {6,12}
    do cd ${f}_${g}_out
    perl ASV_depth_calculation.pl -t ~/path/date_testing/${f}_${g}_out/ASV2Taxonomy/OUT_asvTaxonomyTable.txt -c ~/revamp/outdir/path/dada2/ASVs_counts.tsv -d ~/file/from/btab_date_manipulations.pl/ASV_depthOfTaxonomicInference.txt -o ~/path/originalREVAMPrun/OUT_asvTaxonomyTable.txt > ASV_stats_calculations.txt
    cd ../
    done
  done
```

### Clean up final files

```
for f in {1995..2022}
  do for g in {6,12}
    do cp ${f}_${g}_out/ASV_stats_calculations.txt ASV_stats_calculations_${f}_${g}.txt
    done
  done

for f in {1995..2022}
  do for g in {6,12}
    do echo ${f}_${g} >> ASV_stats_ALLIN.txt
    cat ASV_stats_calculations_${f}_${g}.txt >> ASV_stats_ALLIN.txt
    done
  done
```
