#!/bin/bash

# The list of old file names
old_files=("A1-S13.csv" "A2-S14.csv" "A3-S15.csv" "A4-S16.csv" "A5-S17.csv" \
    "B1-S18.csv" "B2-S19.csv" "B3-S20.csv" "B4-S21.csv" "B5-S22.csv" \
    "C1-S23.csv" "C2-S24.csv" "C3-S25.csv" "C4-S26.csv" "C5-S27.csv" \
    "D1-S28.csv" "D2-S29.csv" "D3-S30.csv" "D4-S31.csv" "D5-S32.csv" \
    "E1-S33.csv" "E2-S34.csv" "E3-S35.csv" "E4-S36.csv" "E5-S37.csv" \
    "F1-S38.csv" "F2-S39.csv" "F3-S40.csv" "F4-S41.csv" "F5-S42.csv" \
    "G1-S43.csv" "G2-S44.csv" "G3-S45.csv" "G4-S46.csv" "G5-S47.csv" \
    "H1-S48.csv" "H2-S49.csv" "H3-S50.csv" "H4-S51.csv" "H5-S52.csv")

# The list of new file names
new_files=('noFMT_1_duodenum' 'noFMT_1_jejunum' 'noFMT_1_ileum' 'noFMT_1_cecum' 'noFMT_1_colon' \
    'noFMT_2_duodenum' 'noFMT_2_jejunum' 'noFMT_2_ileum' 'noFMT_2_cecum' 'noFMT_2_colon' \
    'noFMT_3_duodenum' 'noFMT_3_jejunum' 'noFMT_3_ileum' 'noFMT_3_cecum' 'noFMT_3_colon' \
    'noFMT_4_duodenum' 'noFMT_4_jejunum' 'noFMT_4_ileum' 'noFMT_4_cecum' 'noFMT_4_colon' \
    'donor_1_duodenum' 'donor_1_jejunum' 'donor_1_ileum' 'donor_1_cecum' 'donor_1_colon' \
    'donor_2_duodenum' 'donor_2_jejunum' 'donor_2_ileum' 'donor_2_cecum' 'donor_2_colon' \
    'donor_3_duodenum' 'donor_3_jejunum' 'donor_3_ileum' 'donor_3_cecum' 'donor_3_colon' \
    'donor_4_duodenum' 'donor_4_jejunum' 'donor_4_ileum' 'donor_4_cecum' 'donor_4_colon')

# Iterate over the lists
for i in "${!old_files[@]}"
do
    # Rename the file
    mv "${old_files[i]}" "${new_files[i]}"
done



#!/bin/bash

# The list of files to add the headers to
files=('noFMT_1_duodenum' 'noFMT_1_jejunum' 'noFMT_1_ileum' 'noFMT_1_cecum' 'noFMT_1_colon' \
    'noFMT_2_duodenum' 'noFMT_2_jejunum' 'noFMT_2_ileum' 'noFMT_2_cecum' 'noFMT_2_colon' \
    'noFMT_3_duodenum' 'noFMT_3_jejunum' 'noFMT_3_ileum' 'noFMT_3_cecum' 'noFMT_3_colon' \
    'noFMT_4_duodenum' 'noFMT_4_jejunum' 'noFMT_4_ileum' 'noFMT_4_cecum' 'noFMT_4_colon' \
    'donor_1_duodenum' 'donor_1_jejunum' 'donor_1_ileum' 'donor_1_cecum' 'donor_1_colon' \
    'donor_2_duodenum' 'donor_2_jejunum' 'donor_2_ileum' 'donor_2_cecum' 'donor_2_colon' \
    'donor_3_duodenum' 'donor_3_jejunum' 'donor_3_ileum' 'donor_3_cecum' 'donor_3_colon' \
    'donor_4_duodenum' 'donor_4_jejunum' 'donor_4_ileum' 'donor_4_cecum' 'donor_4_colon')

# Iterate over the list of files
for file in "${files[@]}"
do
    # Add the headers to the first row of the file
    awk -v file_name="$file" 'BEGIN{print "gene," file_name} {print}' $file > temp && mv temp $file
done