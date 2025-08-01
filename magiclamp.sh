#!/bin/bash
eval "$(/home/ark/miniconda3/bin/conda shell.bash hook)"
conda activate base  # Activate the base environment where `boto3` is installed

exec > >(tee -i /home/ark/MAB/magiclamp/magiclamp_looper.log)
exec 2>&1

eval "$(/home/ark/miniconda3/bin/conda shell.bash hook)"
conda activate base  # Activate the base environment where `boto3` is installed

KEY=$1
ID=$KEY
DIR=/home/ark/MAB/magiclamp/${ID}
OUT=/home/ark/MAB/magiclamp/completed/${ID}-results


name=$(grep 'Name' ${DIR}/form-data.txt | cut -d ' ' -f2)
email=$(grep 'Email' ${DIR}/form-data.txt | cut -d ' ' -f2)
option=$(grep 'Option' ${DIR}/form-data.txt | cut -d ' ' -f2)

# Set PATH to include Conda and script locations
export PATH="/home/ark/miniconda3/bin:/usr/local/bin:/usr/bin:/bin:/home/ark/MAB/bin/magiclamp-local:$PATH"
eval "$(/home/ark/miniconda3/bin/conda shell.bash hook)"
conda activate magiclamp

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate Conda environment."
    exit 1
fi
sleep 5

# Rename files
for file in ${DIR}/*.f*; do
    [[ ${file} == ${DIR}/form-data.txt ]] || mv ${file} ${file%.*}.fa
    /home/ark/MAB/bin/BagOfTricks/header-format.py -file ${file%.*}.fa -out ${file%.*}.fxa -char '|' -rep '-'
done

if [[ ${option} == "Custom" ]]; then

    mkdir ${DIR}/HMMs
    mv ${DIR}/*.hmm ${DIR}/HMMs/
    mv ${DIR}/*.HMM ${DIR}/HMMs/

    for file in ${DIR}/HMMs/*; do
      mv ${file} ${file%.*}.hmm
    done

    echo /home/ark/bin/MAB/MagicLamp/MagicLamp.py HmmGenie -bin_dir ${DIR} -bin_ext fxa -out ${OUT} -t 4 -hmm_dir ${DIR}/HMMs -hmm_ext hmm
    /home/ark/MAB/bin/MagicLamp/MagicLamp.py HmmGenie -bin_dir ${DIR} -bin_ext fxa -out ${OUT} -t 4 -hmm_dir ${DIR}/HMMs -hmm_ext hmm
elif [[ ${option} == "FeGenie" ]]; then
    echo /home/ark/MAB/bin/MagicLamp/MagicLamp.py FeGenie -bin_dir ${DIR} -bin_ext fxa -out ${OUT} -t 4 -hmm_dir ${DIR}/HMMs -hmm_ext hmm
    /home/ark/MAB/bin/MagicLamp/MagicLamp.py FeGenie -bin_dir ${DIR} -bin_ext fxa -out ${OUT} -t 4
elif [[ ${option} == "LithoGenie" ]]; then
    echo /home/ark/MAB/bin/MagicLamp/MagicLamp.py LithoGenie -bin_dir ${DIR} -bin_ext fxa -out ${OUT} -t 4 -hmm_dir ${DIR}/HMMs -hmm_ext hmm
    /home/ark/MAB/bin/MagicLamp/MagicLamp.py LithoGenie -bin_dir ${DIR} -bin_ext fxa -out ${OUT} -t 4
else
    echo /home/ark/MAB/bin/MagicLamp/MagicLamp.py OmniGenie -bin_dir ${DIR} -bin_ext fxa -out ${OUT} -t 4 -genie ${option}
    /home/ark/MAB/bin/MagicLamp/MagicLamp.py OmniGenie -bin_dir ${DIR} -bin_ext fxa -out ${OUT} -t 4 -genie ${option}
fi

# **************************************************************************************************
# **************************************************************************************************
# **************************************************************************************************
if [ $? -ne 0 ]; then
    echo "Error: MagicLamp failed."
    conda deactivate
    exit 1
fi
conda deactivate
sleep 5

# Archive results
mv /home/ark/MAB/magiclamp/completed/${ID}-results ./${ID}-results
tar -cf ${ID}-results.tar ${ID}-results && gzip ${ID}-results.tar

# Upload results to S3 and generate presigned URL
results_tar="${ID}-results.tar.gz"
s3_key="${ID}-results.tar.gz"
python3 /home/ark/MAB/bin/magiclamp-local/push.py --bucket binfo-dump --output_key ${s3_key} --source ${results_tar}
url=$(python3 /home/ark/MAB/bin/magiclamp-local/gen_presign_url.py --bucket binfo-dump --key ${s3_key} --expiration 86400)

mv ${ID}-results.tar.gz /home/ark/MAB/magiclamp/completed/${ID}-results.tar.gz
rm -rf ./${ID}-results


# Send email
python3 /home/ark/MAB/bin/magiclamp-local/send_email.py \
    --sender ark@midauthorbio.com \
    --recipient ${email} \
    --subject "Your MagicLamp Results!" \
    --body "Hi ${name},

    Your MagicLamp results are available for download using the link below. The link will expire in 24 hours.

    ${url}

    Please visit https://github.com/Arkadiy-Garber/MagicLamp for documentation.

    Please reach out to ark@midauthorbio.com, or send us a note on https://midauthorbio.com/#contact if you have any questions.

    Thanks!
    MAB Team"

if [ $? -ne 0 ]; then
    echo "Error: send_email.py failed."
    conda deactivate
    exit 1
fi

sleep 5

#sudo rm -rf ${DIR}

conda deactivate
echo "MagicLamp completed successfully."



