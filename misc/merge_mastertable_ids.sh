panel=Maanasa_americas_reheaded_filtered

# Save samples' IDs to a file
bcftools view -h $panel.vcf | sed 's/\t/\n/g' |grep -v '#' \
 |tail -n+9 | cut -f1 -d '_' >names.$panel.txt

echo ${selected[@]} |sed 's/ /\n/g' > names.$name.txt

echo "ID" > names.$name.$panel.txt
cat names.$name.txt names.$panel.txt >>names.$name.$panel.txt

# Merge a 'master table' to the IDs
