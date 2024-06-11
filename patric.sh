  GNU nano 4.8                                                                                                      script.sh                                                                                                                #!/usr/bin/bash
############### Scripts: get all list of genome PATRIC
## Autor: Leandro de Mattos
### Step1:

#for i in `cat Genomes_id`;
#do mkdir "$i"
#done

## Step2
for i in `cat Genomes_id`;
#do cd "$i"
do
wget -qN "ftp://ftp.patricbrc.org/genomes/$i/$i*.gb";

done
#wget -qN "ftp://ftp.patricbrc.org/genomes/$i/$i.*";
