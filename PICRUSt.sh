## run the first part of the lightPollution_picrust.R first

export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate picrust2

#https://github.com/picrust/picrust2/wiki/PICRUSt2-Tutorial-(v2.3.0-beta)#place-reads-into-reference-tree

place_seqs.py -s light_pollution_ASV.fasta -o light_pollution.tre -p 64 --intermediate intermediate/place_seqs 

#Warning - 75 input sequences aligned poorly to reference sequences (--min_align option specified a minimum proportion of 0.8 aligning to reference sequences). These input sequences will not be placed and will be excluded from downstream steps.
#This is the set of poorly aligned input sequences to be excluded: ASV18860, ASV13779, ASV18182, ASV18290, ASV18408, ASV18838, ASV18469, ASV18544, ASV17984, ASV11742, ASV18168, ASV18658, ASV4806, ASV18640, ASV4392, ASV16589, ASV18807, ASV18465, ASV18782, ASV16565, ASV18749, ASV18143, ASV17873, ASV18468, ASV10597, ASV15053, ASV18855, ASV18144, ASV17756, ASV17763, ASV18696, ASV18170, ASV5179, ASV18784, ASV18033, ASV17110, ASV15378, ASV17764, ASV18055, ASV18101, ASV18642, ASV18596, ASV17435, ASV14055, ASV16851, ASV13908, ASV18646, ASV10579, ASV18427, ASV18454, ASV18190, ASV18318, ASV9748, ASV3107, ASV1182, ASV18464, ASV18396, ASV16848, ASV18301, ASV18645, ASV15957, ASV8863, ASV12685, ASV18806, ASV18815, ASV16663, ASV7381, ASV6178, ASV13778, ASV8365, ASV15559, ASV18460, ASV16342, ASV18227, ASV18399

hsp.py -i 16S -t light_pollution.tre -o marker_predicted_and_nsti.tsv.gz -p 64 -n
hsp.py -i KO -t light_pollution.tre -o KO_predicted.tsv.gz -p 64
metagenome_pipeline.py -i light_pollution_ASV_picrust -m marker_predicted_and_nsti.tsv.gz -f KO_predicted.tsv.gz  -o KO_metagenome_out --strat_out

## run the second part of the lightPollution_picrust.R to analyze the results.
