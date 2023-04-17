#!/bin/bash
### Automating smash on AWS: use representative samples from each individual to check for swaps across batches
### micah 2022.02.07

# exit on error
set -o errexit

## 0) Use "rep_samples_s3_paths.txt" to specify which batches to run, their IDs and their location on s3
project_dir='/home/ec2-user/micah/2023/stanford_smashqc'
s3_paths=$(cut -f2 -d " " ${project_dir}/processing/scripts/rep_samples_s3_paths_20230214.txt)

echo 'rep_samples_s3_paths_20230214.txt:'
cat ${project_dir}/processing/scripts/rep_samples_s3_paths_20230214.txt

# turn space sep string into comma sep string
function filelistjoin { local IFS="$1"; shift; echo "$*"; };

# load smash docker image
zcat /home/ec2-user/docker/Images/smash.tar.gz | docker load

#n_batches=$(echo $batches | wc -w)

# 1) Sync the fastq data from s3 (will need to rm these when run is complete to make space for more)
batch_id='rep_samps_ALL'
echo $batch_id

mkdir -pv ${project_dir}/data/${batch_id}
cd ${project_dir}/data/${batch_id}

for s3_path in $s3_paths; 
do
	echo $s3_path
	echo `dirname $s3_path`
	echo `basename $s3_path`
	aws s3 sync `dirname $s3_path` . --exclude "*" --include `basename $s3_path`
done

# 2) Create subdir structure
mkdir -pv ${project_dir}/processing/${batch_id}/star
mkdir -pv ${project_dir}/processing/${batch_id}/smash/bams
mkdir -pv ${project_dir}/processing/${batch_id}/smash/output

# 3) Align each fastq in batch using star
cd ${project_dir}/processing/${batch_id}

for r1p in ${project_dir}/data/${batch_id}/*L00?_R1*.fastq.gz; # for each R1 fastq
 do 
 	r1f=`basename ${r1p}`;
         sid=${r1f%%_S*_L00?*};
 	if [ ! -e star/${sid}/${sid}_Log.std.out ];
         then 
 		echo ${sid};
                 mkdir -pv star/${sid}/;
                 r1list=(${project_dir}/data/${batch_id}/${sid}_*R1*.fastq.gz);
                 r1files=`filelistjoin , ${r1list[@]}`;
                 r2list=(${project_dir}/data/${batch_id}/${sid}_*R2*.fastq.gz);
                 r2files=`filelistjoin , ${r2list[@]}`;
                 echo ${r1files} ${r2files};
                 /home/ec2-user/tools/STAR/STAR-2.7.9a/bin/Linux_x86_64_static/STAR \
 			--genomeDir /home/ec2-user/References/homo_sapiens/GRCh38/ensembl_107/star_index_2.7.9a \
                         --genomeLoad LoadAndKeep \
                         --readFilesIn ${r1files} ${r2files} \
 			--readFilesCommand zcat \
                         --runThreadN 64 \
                         --outStd BAM_Unsorted \
                         --outSAMtype BAM Unsorted \
                         --outSAMmode Full \
                         --outSAMattributes All \
                         --outSAMunmapped Within \
                         --outTmpDir /home/ec2-user/micah/2023/stanford_smashqc/processing/tmp/starTmp \
                         --outFileNamePrefix ${project_dir}/processing/${batch_id}/star/${sid}/${sid}_ \
                         --outFilterMultimapNmax 999 \
                         --outSAMprimaryFlag AllBestScore \
                         --quantMode GeneCounts > ${project_dir}/processing/${batch_id}/star/${sid}/${sid}_GRCh38-107_S279a_unsorted.bam \
                         --readMapNumber 5000000;
 	fi;
 done

# 4) Sort and index bams
cd ${project_dir}/processing/${batch_id}
for ubp in star/*/*.bam; do \
        ubdir=`dirname ${ubp}`; \
        ubf=`basename ${ubp}`; \
        ubpre=${ubf%%_unsorted.bam}; \
        echo ${ubp} ${ubdir} ${ubf} ${ubpre}\n; \
        /home/ec2-user/tools/samtools/samtools-1.10/samtools sort -@ 64 -m 2G -o ${ubdir}/${ubpre}_sorted.bam ${ubp}; \
done

for sbf in star/*/*_sorted.bam; do \
        echo ${sbf}; \
        /home/ec2-user/tools/samtools/samtools-1.10/samtools index -@ 64 ${sbf}; \
done

# 5) Link star output sorted bams to smash bams dir
cd ${project_dir}/processing/${batch_id}/smash/bams

ln -v ${project_dir}/processing/${batch_id}/star/*/*_sorted.bam .
ln -v ${project_dir}/processing/${batch_id}/star/*/*_sorted.bam.bai .

# 6) Run smash
docker run \
        --rm -it \
        -v /home/ec2-user/tools/SMaSH-master:/home/ec2-user/tools/SMaSH-master \
        -v ${project_dir}/processing/${batch_id}/smash/bams:/smash/bams \
	-v ${project_dir}/processing/${batch_id}/smash/output/:/smash/output \
        -w /smash/bams \
        ff13285f3b4f \
        /home/ec2-user/tools/SMaSH-master/SMaSH.py \
                -i /home/ec2-user/tools/SMaSH-master/snps_GRCh38.vcf \
                -bam *.bam \
                -output_dir ../output/

# 7) Check if complete, if yes then rm all tmp and sequence files and continue to next batch, if not then break loop
if [ -e ${project_dir}/processing/${batch_id}/smash/output/pval_out.txt ]
then
	rm -vrf ${project_dir}/data/${batch_id}/*.fastq* || true
	rm -vrf ${project_dir}/processing/${batch_id}/star/*/*.ba? || true
	rm -vrf ${project_dir}/processing/${batch_id}/smash/bams/*
	rm -vrf ${project_dir}/processing/tmp/starTmp/* || true
else
	echo "SMaSH output is missing!"
	break
fi
