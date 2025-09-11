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

# --- helper: detect if FASTA is nucleotide (ACGTN only in first non-header seq chunk)
is_nucleotide_fasta () {
  local fasta="$1"
  # grab first non-header line(s), strip whitespace, limit to 2000 chars, then test
  local sample
  sample=$(grep -v '^>' "$fasta" | tr -d ' \t\r\n' | head -c 2000 || true)
  if [[ -z "$sample" ]]; then
    # empty sequence? assume nucleotide to be conservative
    return 0
  fi
  if [[ "$sample" =~ ^[ACGTNacgtn]+$ ]]; then
    return 0  # nucleotide
  else
    return 1  # protein
  fi
}

name=$(grep 'Name' ${DIR}/form-data.txt | cut -d ' ' -f2)
email=$(grep 'Email' ${DIR}/form-data.txt | cut -d ' ' -f2)
# Multiple references supported: collect all third fields on lines beginning with "Input"
mapfile -t REF_FILES < <(awk '/^Ref/ {print $3}' "${DIR}/form-data.txt")

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

echo /home/ark/MAB/bin/crystalomics-local/cif-peptide-extract.py -cif ${DIR}/${cif} -faa ${OUT}/${cif%.*}.faa -txt ${OUT}/${cif%.*}.txt
/home/ark/MAB/bin/crystalomics-local/cif-peptide-extract.py -cif ${DIR}/${cif} -faa ${OUT}/${cif%.*}.faa -txt ${OUT}/${cif%.*}.txt

# --- loop over each reference in form-data.txt
for ref_rel in "${REF_FILES[@]}"; do
  ref_path="${DIR}/${ref_rel}"
  ref_base="$(basename "${ref_rel}")"
  ref_label="${ref_base%.*}"
  out_blast="${OUT}/${ref_label}.blast"

  echo "Preparing DB for reference: ${ref_rel}"

  if is_nucleotide_fasta "${ref_path}"; then
    echo "Detected nucleotide DB → makeblastdb -dbtype nucl + tblastn"
    makeblastdb -in "${ref_path}" -dbtype nucl -out "${ref_path}" 1>&2
    tblastn -num_threads 16 \
            -out "${out_blast}" -outfmt 6 \
            -query "${OUT}/${cif%.*}.faa" \
            -db "${ref_path}" \
            -max_target_seqs 1
  else
    echo "Detected protein DB → makeblastdb -dbtype prot + blastp"
    makeblastdb -in "${ref_path}" -dbtype prot -out "${ref_path}" 1>&2
    blastp -num_threads 16 \
           -out "${out_blast}" -outfmt 6 \
           -query "${OUT}/${cif%.*}.faa" \
           -db "${ref_path}" \
           -max_target_seqs 1
  fi

  BLAST_FILES+=("${out_blast}")
done

/home/ark/MAB/bin/crystalomics-local/blast2summary.v2.py -f ${OUT}/${cif%.*}.faa -o ${OUT}/${cif%.*}.ref.summary.csv -b "${BLAST_FILES[@]}"


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

