## Version - inStrain/1.5.3
## Script for inStrain strain profiling and comparison

## Strain profiling module

for i in $(cat samplenames.txt)
do

    ## bowtie2/2.3.4
    bowtie2 -p 30 -x InStrain_DB/reps.fasta.bt2 -1 sample_R1.fastq -2 sample_R2.fastq > "$i".sam
    

    ## samtools/1.10
    samtools view "$i".sam -b -o "$i".bam -@ 19
    rm "$i".sam
	
    ## inStrain Profiling 
   
    inStrain profile "$i".bam InStrain_DB/reps.fasta -o "$i"_instrainout -p 30 -g InStrain_DB/Genes/reps.genes.fna -s InStrain_DB/genomes.stb --database_mode
 

    ## Removing empty instrainout folders
    [ "$(ls -A "$i"_instrainout/output/)" ] && echo "Not Empty" || rm -r "$i"_instrainout/

done

## inStrain profile comparison

inStrain compare -ani 0.99999 -cov 0.5 -i *_instrainout/ -s instrain_built_database.stb -p 25 -o output_dir --database_mode


