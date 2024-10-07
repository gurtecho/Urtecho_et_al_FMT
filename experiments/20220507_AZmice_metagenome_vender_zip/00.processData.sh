for eachS in `cat assembly.list`
do
	IFS='_'; arrIN=($eachS); unset IFS;
	supplier=$arrIN
	echo $supplier
#	./bin/processAnnotation_COG.py ./output_12_reAssemblyFunctionalAnnotation/$supplier\_COG/$eachS\.COG.txt \
#					./processed_data/$eachS\.COG.txt 0.0001
#	./bin/processAnnotation_CAZ.py ./output_12_reAssemblyFunctionalAnnotation/$supplier\_CAZy/$eachS\.CAZy.txt \
#					./processed_data/$eachS\.CAZy_GH.txt GH 0.0001
#	./bin/processAnnotation_KEGG.py ./output_12_reAssemblyFunctionalAnnotation/$supplier\_KEGG/$eachS\.KEGG.txt \
#					./processed_data/$eachS\.KEGG.txt 0.0001
#	bin/processAnnotation_SignalP.py ./output_12_reAssemblyFunctionalAnnotation/$supplier\_SignalIP/$eachS\_summary.signalp5 \
#					./processed_data/$eachS\.signalP.txt $eachS
done
#cat ./processed_data/*.COG.txt > assembly.COG.txt
#cat ./processed_data/*.CAZy_GH.txt > assembly.CAZy_GH.txt
#cat ./processed_data/*.KEGG.txt > assembly.KEGG.txt
#cat ./processed_data/*.signalP.txt > assembly.signalP.txt

