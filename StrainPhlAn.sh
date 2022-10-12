## Script for StrainPhlAn - identify strains of a species

## Version -  metaphlan2/3.0

## consensus-marker files -  input for StrainPhlAn

for i in $(cat samplenames.txt)
do
    sample2markers.py -i MetaPhlanAn3_Outputs/"$i"_metaphlan3.sam.bz2 -o ConsensusMarkers --nproc 10
done

## extracting clade names

strainphlan -s ConsensusMarkers/*.pkl -d database/mpa_v30_CHOCOPhlAn_201901.pkl --print_clades_only -o . --nproc 10 > clades.txt

## Extract the markers from MetaPhlAn database
## Build the multiple sequence alignment and the phylogenetic tree

for i in $(grep 's__' clades.txt | sed 's/.*s__/s__/' | sed 's/: .*//')
do
    extract_markers.py -c "$i" -d database/mpa_v30_CHOCOPhlAn_201901.pkl -o CladeMarkers/
    strainphlan -s ConsensusMarkers/*.pkl -m CladeMarkers/"$i".fna -o Output/ -c "$i" -d database/mpa_v30_CHOCOPhlAn_201901.pkl --nproc 10
done

