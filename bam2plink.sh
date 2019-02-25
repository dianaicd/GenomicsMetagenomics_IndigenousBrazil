panel=Maanasa_mask1_flip


while read botocudo
    do
    echo $botocudo
    ind=$(basename $botocudo .bam)
    bam2plink.py bamfile=$botocudo plinkpref=$panel MinBQ=20 indname=$botocudo popname=Botocudo doCns=F trim=5 MinMQ=30 &
done < botocudos.txt

{ sleep 5; echo waking up after 5 seconds; } &
{ sleep 1; echo waking up after 1 second; } &
  wait
  echo all jobs are done!

while read botocudo
    do
    echo $botocudo
    ind=$(basename $botocudo .bam)
    bam2plink.py bamfile=$botocudo plinkpref=$panel MinBQ=20 indname=$botocudo popname=$botocudo doCns=F trim=5 MinMQ=30 &
done < fromvictor.txt

{ sleep 5; echo waking up after 5 seconds; } &
{ sleep 1; echo waking up after 1 second; } &
  wait
  echo all jobs are done!

# Merge plink binary files
 ln -s ~/archive/Panels/fromVictor/Maanasa_mask1_flip.{bed,bim,fam} ./
maanasa=Maanasa_mask1_flip

panel=Maanasa_mask1_flip

while read line 
do
    boto=$(basename $line .bam)
    echo "-------------------- Merging $boto to $panel ---------------------"

    plink --bfile $panel --bmerge plink_botocudo/${boto}.bam_Botocudo_${maanasa} --make-bed --out ${maanasa}_${boto} --allow-no-sex
    rm $panel.*
    panel=${maanasa}_${boto}
done < botocudos.txt

 ln -s ~/archive/Panels/fromVictor/Maanasa_mask1_flip.{bed,bim,fam} ./
maanasa=Maanasa_mask1_flip

panel=Maanasa_mask1_flip

while read line 
do
    boto=$(basename $line .bam)
    echo "-------------------- Merging $boto to $panel ---------------------"

    plink --bfile $panel --bmerge plink_botocudo/${boto}.bam_Botocudo_${maanasa} --make-bed --out ${maanasa}_${boto} --allow-no-sex
    rm $panel.*
    panel=${maanasa}_${boto}
done < botocudos.txt

 ln -s ~/archive/Panels/fromVictor/Maanasa_mask1_flip.{bed,bim,fam} ./
maanasa=Maanasa_mask1_flip

while read line 
do
    boto=$(basename $line .bam)
    echo "-------------------- Merging $boto to $panel ---------------------"

    plink --bfile $panel --bmerge fromVictor/${boto}.bam_${boto}.bam_${maanasa} --make-bed --out ${maanasa}_${boto} --allow-no-sex
    rm $panel.*
    panel=${maanasa}_${boto}
done < fromvictor.txt
