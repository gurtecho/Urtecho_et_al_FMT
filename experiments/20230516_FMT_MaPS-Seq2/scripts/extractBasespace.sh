for i in *ds.* 
do	
	cd $i
	mv *.fastq.gz ..
	cd ..
done

rm -r *ds.*

mv * ..
cd ..

for file in *_S*R1_001.fastq.gz; do mv "$file" "${file/_S*_001.fastq.gz/}_R1.fastq.gz"; done
for file in *_S*R2_001.fastq.gz; do mv "$file" "${file/_S*_001.fastq.gz/}_R2.fastq.gz"; done

