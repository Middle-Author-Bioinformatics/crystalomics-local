#!/bin/bash
eval "$(/home/ark/miniconda3/bin/conda shell.bash hook)"
conda activate base  # Activate the base environment where `boto3` is installed

exec > >(tee -i /home/ark/MAB/crystalomics/crystalomics_looper.log)
exec 2>&1

eval "$(/home/ark/miniconda3/bin/conda shell.bash hook)"
conda activate base  # Activate the base environment where `boto3` is installed

KEY=$1
ID=$KEY
DIR=/home/ark/MAB/crystalomics/${ID}
OUT=/home/ark/MAB/crystalomics/completed/${ID}-results
mkdir -p ${OUT}


name=$(grep 'Name' ${DIR}/form-data.txt | cut -d ' ' -f2)
email=$(grep 'Email' ${DIR}/form-data.txt | cut -d ' ' -f2)
reference=$(grep 'Input' ${DIR}/form-data.txt | cut -d ' ' -f3)
cif=$(grep 'CIF' ${DIR}/form-data.txt | cut -d ' ' -f3)

# Set PATH to include Conda and script locations
export PATH="/home/ark/miniconda3/bin:/usr/local/bin:/usr/bin:/bin:/home/ark/MAB/bin/crystalomics-local:$PATH"
#eval "$(/home/ark/miniconda3/bin/conda shell.bash hook)"
#conda activate crystalomics

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate Conda environment."
    exit 1
fi
sleep 5

# **************************************************************************************************
# **************************************************************************************************
# **************************************************************************************************

for i in ${cif}; do
    echo "Processing CIF file: ${i}"
    /home/ark/MAB/bin/crystalomics-local/cif-peptide-extract.py -cif ${DIR}/${i} -faa ${OUT}/${i%.*}.faa -txt ${OUT}/${i%.*}.txt
    makeblastdb -dbtype prot -in ${DIR}/${reference} -out ${DIR}/${reference}
    blastp -num_threads 16 -out ${OUT}/${i%.*}.ref.blast -outfmt 6 -query ${OUT}/${i%.*}.faa -db ${DIR}/${reference} -max_target_seqs 1
    /home/ark/MAB/bin/crystalomics-local/blast2summary.py -db ${DIR}/${reference} -b ${OUT}/${i%.*}.ref.blast -f ${OUT}/${i%.*}.faa -o ${OUT}/${i%.*}.ref.summary.csv
    diamond blastp --threads 16 -d /home/ark/databases/nr.dmnd -q ${OUT}/${i%.*}.faa -o ${OUT}/${i%.*}.nr.blastp -f 6 qseqid sseqid pident length evalue bitscore stitle qseq sseq --max-target-seqs 10 --evalue 10
    /home/ark/MAB/bin/crystalomics-local/nr2summary.py -b1 ${OUT}/${i%.*}.nr.blastp -b2 ${OUT}/${i%.*}.ref.blast -o ${OUT}/${i%.*}.nr.summary.tsv
done


# **************************************************************************************************
# **************************************************************************************************
# **************************************************************************************************
if [ $? -ne 0 ]; then
    echo "Error: Crystalomics failed."
#    conda deactivate
    exit 1
fi
#conda deactivate
#sleep 5

# Archive results
mv /home/ark/MAB/crystalomics/completed/${ID}-results ./${ID}-results
tar -cf ${ID}-results.tar ${ID}-results && gzip ${ID}-results.tar

# Upload results to S3 and generate presigned URL
results_tar="${ID}-results.tar.gz"
s3_key="${ID}-results.tar.gz"
python3 /home/ark/MAB/bin/crystalomics-local/push.py --bucket binfo-dump --output_key ${s3_key} --source ${results_tar}
url=$(python3 /home/ark/MAB/bin/crystalomics-local/gen_presign_url.py --bucket binfo-dump --key ${s3_key} --expiration 86400)

mv ${ID}-results.tar.gz /home/ark/MAB/crystalomics/completed/${ID}-results.tar.gz
rm -rf ./${ID}-results


# Send email
python3 /home/ark/MAB/bin/crystalomics-local/send_email.py \
    --sender ark@midauthorbio.com \
    --recipient ${email} \
    --subject "Your CIF peptides!" \
    --body "Hi ${name},

    Your CIF peptide results are available for download using the link below. The link will expire in 24 hours.

    ${url}

    Please reach out to agarber4@asu.com if you have any questions.

    Thanks!
    Your friendly neighborhood bioinformatician 🕸️"

if [ $? -ne 0 ]; then
    echo "Error: send_email.py failed."
#    conda deactivate
    exit 1
fi

sleep 5

#sudo rm -rf ${DIR}

#conda deactivate
echo "Crystalomics completed successfully."



