for i in ds.*;

do 

cd $i

wc -l *_R1*

cd ..


done >> ./sampleCounts.txt

