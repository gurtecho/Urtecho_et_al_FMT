#!/bin/bash

# The list of old file names
old_files=("postFMTs1scecum_S27.csv" "postFMTs1scolon_S28.csv" "postFMTs1sduodenum_S24.csv" "postFMTs1sileum_S26.csv" "postFMTs1sjejunum_S25.csv" \
    "postFMTs2scecum_S32.csv" "postFMTs2scolon_S33.csv" "postFMTs2sduodenum_S29.csv" "postFMTs2sileum_S31.csv" "postFMTs2sjejunum_S30.csv" \
    "postFMTs3scecum_S37.csv" "postFMTs3scolon_S38.csv" "postFMTs3sduodenum_S34.csv" "postFMTs3sileum_S36.csv" "postFMTs3sjejunum_S35.csv" \
    "postFMTs4scecum_S42.csv" "postFMTs4scolon_S43.csv" "postFMTs4sduodenum_S39.csv" "postFMTs4sileum_S41.csv" "postFMTs4sjejunum_S40.csv")

# The list of new file names
new_files=('postFMT_1_cecum' 'postFMT_1_colon' 'postFMT_1_duodenum' 'postFMT_1_ileum' 'postFMT_1_jejunum' \
    'postFMT_2_cecum' 'postFMT_2_colon' 'postFMT_2_duodenum' 'postFMT_2_ileum' 'postFMT_2_jejunum' \
    'postFMT_3_cecum' 'postFMT_3_colon' 'postFMT_3_duodenum' 'postFMT_3_ileum' 'postFMT_3_jejunum' \
    'postFMT_4_cecum' 'postFMT_4_colon' 'postFMT_4_duodenum' 'postFMT_4_ileum' 'postFMT_4_jejunum')

# Iterate over the lists
for i in "${!old_files[@]}"
do
    # Rename the file
    mv "${old_files[i]}" "${new_files[i]}"
done



#!/bin/bash

# The list of files to add the headers to
files=('postFMT_1_cecum' 'postFMT_1_colon' 'postFMT_1_duodenum' 'postFMT_1_ileum' 'postFMT_1_jejunum' \
    'postFMT_2_cecum' 'postFMT_2_colon' 'postFMT_2_duodenum' 'postFMT_2_ileum' 'postFMT_2_jejunum' \
    'postFMT_3_cecum' 'postFMT_3_colon' 'postFMT_3_duodenum' 'postFMT_3_ileum' 'postFMT_3_jejunum' \
    'postFMT_4_cecum' 'postFMT_4_colon' 'postFMT_4_duodenum' 'postFMT_4_ileum' 'postFMT_4_jejunum')

# Iterate over the list of files
for file in "${files[@]}"
do
    # Add the headers to the first row of the file
    awk -v file_name="$file" 'BEGIN{print "gene," file_name} {print}' $file > temp && mv temp $file
done