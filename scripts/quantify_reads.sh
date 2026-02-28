#!/bin/bash

for i in ../reads/SRR105516{57..65};
do

read=`basename ${i}`

salmon quant -i transcripts_index \
	-l A \
	-r ../reads/${read}.fastq \
	-p 6 \
	--validateMappings \
	--gcBias \ 
	--seqBias \
	-o quants/${read} 

done

for i in quants/SRR105516{57..65};

do

echo "Mapping rate of ${i} is: $(grep -oP "Mapping rate = \K[\d.]+" ${i}/logs/salmon_quant.log)"

done
