#!/bin/bash

# Usage help
usage() {
  echo "Usage: $0 <input_json_file> --data-type <joint|exome|genome> -o <output_directory>"
  exit 1
}

# Parse arguments
input_file=""
data_type=""
output_dir=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--output)
      output_dir="$2"
      shift 2
      ;;
    --data-type)
      data_type="$2"
      shift 2
      ;;
    -*)
      echo "Unknown option: $1"
      usage
      ;;
    *)
      input_file="$1"
      shift
      ;;
  esac
done

# Validate input
if [[ -z "$input_file" || -z "$data_type" || -z "$output_dir" ]]; then
  usage
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Exit if variant not found
if ! jq -e '.data.variant != null' "$input_file" > /dev/null; then
  echo "⚠️  Skipping (variant not found): $input_file"
  exit 0
fi

# Create output file path
base_name=$(basename "$input_file" .json)
output_file="$output_dir/${base_name}.tsv"

# Header fields
header_fields='["variant_id", "gene_symbol", "pos", "rsids", "transcript_id", "transcript_version", "hgvs", "hgvsp", "consequence", 
  "'"$data_type"'_ac", "'"$data_type"'_ac_hemizygote_count", "'"$data_type"'_ac_homozygote_count", "'"$data_type"'_an",
  "'"$data_type"'_afr_ac", "'"$data_type"'_afr_an", "'"$data_type"'_afr_ac_hemizygote_count", "'"$data_type"'_afr_ac_homozygote_count",
  "'"$data_type"'_amr_ac", "'"$data_type"'_amr_an", "'"$data_type"'_amr_ac_hemizygote_count", "'"$data_type"'_amr_ac_homozygote_count",
  "'"$data_type"'_asj_ac", "'"$data_type"'_asj_an", "'"$data_type"'_asj_ac_hemizygote_count", "'"$data_type"'_asj_ac_homozygote_count",
  "'"$data_type"'_eas_ac", "'"$data_type"'_eas_an", "'"$data_type"'_eas_ac_hemizygote_count", "'"$data_type"'_eas_ac_homozygote_count",
  "'"$data_type"'_nfe_ac", "'"$data_type"'_nfe_an", "'"$data_type"'_nfe_ac_hemizygote_count", "'"$data_type"'_nfe_ac_homozygote_count",
  "'"$data_type"'_fin_ac", "'"$data_type"'_fin_an", "'"$data_type"'_fin_ac_hemizygote_count", "'"$data_type"'_fin_ac_homozygote_count",
  "'"$data_type"'_mid_ac", "'"$data_type"'_mid_an", "'"$data_type"'_mid_ac_hemizygote_count", "'"$data_type"'_mid_ac_homozygote_count",
  "'"$data_type"'_sas_ac", "'"$data_type"'_sas_an", "'"$data_type"'_sas_ac_hemizygote_count", "'"$data_type"'_sas_ac_homozygote_count",
  "'"$data_type"'_ami_ac", "'"$data_type"'_ami_an", "'"$data_type"'_ami_ac_hemizygote_count", "'"$data_type"'_ami_ac_homozygote_count",
  "'"$data_type"'_remaining_ac", "'"$data_type"'_remaining_an", "'"$data_type"'_remaining_ac_hemizygote_count", "'"$data_type"'_remaining_ac_homozygote_count"
]'
echo "$header_fields" | jq -r '. | @tsv' > "$output_file"

# Write variant row
jq -r --arg dt "$data_type" '
.data.variant as $v |
[
  $v.variant_id,
  ($v.transcript_consequences[0].gene_symbol // "NA"),
  $v.pos,
  ($v.rsids | join(",") // "NA"),
  ($v.transcript_consequences[0].transcript_id // "NA"),
  ($v.transcript_consequences[0].transcript_version // "NA"),
  ($v.transcript_consequences[0].hgvsc // "NA"),
  ($v.transcript_consequences[0].hgvsp // "NA"),
  ($v.transcript_consequences[0].major_consequence // "NA"),

  ($v[$dt].ac // 0),
  ($v[$dt].hemizygote_count // 0),
  ($v[$dt].homozygote_count // 0),
  ($v[$dt].an // 0),

  ($v[$dt].populations[]? | select(.id=="afr") | .ac // 0),
  ($v[$dt].populations[]? | select(.id=="afr") | .an // 0),
  ($v[$dt].populations[]? | select(.id=="afr") | .hemizygote_count // 0),
  ($v[$dt].populations[]? | select(.id=="afr") | .homozygote_count // 0),

  ($v[$dt].populations[]? | select(.id=="amr") | .ac // 0),
  ($v[$dt].populations[]? | select(.id=="amr") | .an // 0),
  ($v[$dt].populations[]? | select(.id=="amr") | .hemizygote_count // 0),
  ($v[$dt].populations[]? | select(.id=="amr") | .homozygote_count // 0),

  ($v[$dt].populations[]? | select(.id=="asj") | .ac // 0),
  ($v[$dt].populations[]? | select(.id=="asj") | .an // 0),
  ($v[$dt].populations[]? | select(.id=="asj") | .hemizygote_count // 0),
  ($v[$dt].populations[]? | select(.id=="asj") | .homozygote_count // 0),

  ($v[$dt].populations[]? | select(.id=="eas") | .ac // 0),
  ($v[$dt].populations[]? | select(.id=="eas") | .an // 0),
  ($v[$dt].populations[]? | select(.id=="eas") | .hemizygote_count // 0),
  ($v[$dt].populations[]? | select(.id=="eas") | .homozygote_count // 0),

  ($v[$dt].populations[]? | select(.id=="nfe") | .ac // 0),
  ($v[$dt].populations[]? | select(.id=="nfe") | .an // 0),
  ($v[$dt].populations[]? | select(.id=="nfe") | .hemizygote_count // 0),
  ($v[$dt].populations[]? | select(.id=="nfe") | .homozygote_count // 0),

  ($v[$dt].populations[]? | select(.id=="fin") | .ac // 0),
  ($v[$dt].populations[]? | select(.id=="fin") | .an // 0),
  ($v[$dt].populations[]? | select(.id=="fin") | .hemizygote_count // 0),
  ($v[$dt].populations[]? | select(.id=="fin") | .homozygote_count // 0),

  ($v[$dt].populations[]? | select(.id=="mid") | .ac // 0),
  ($v[$dt].populations[]? | select(.id=="mid") | .an // 0),
  ($v[$dt].populations[]? | select(.id=="mid") | .hemizygote_count // 0),
  ($v[$dt].populations[]? | select(.id=="mid") | .homozygote_count // 0),

  ($v[$dt].populations[]? | select(.id=="sas") | .ac // 0),
  ($v[$dt].populations[]? | select(.id=="sas") | .an // 0),
  ($v[$dt].populations[]? | select(.id=="sas") | .hemizygote_count // 0),
  ($v[$dt].populations[]? | select(.id=="sas") | .homozygote_count // 0),

  ($v[$dt].populations[]? | select(.id=="ami") | .ac // 0),
  ($v[$dt].populations[]? | select(.id=="ami") | .an // 0),
  ($v[$dt].populations[]? | select(.id=="ami") | .hemizygote_count // 0),
  ($v[$dt].populations[]? | select(.id=="ami") | .homozygote_count // 0),

  ($v[$dt].populations[]? | select(.id=="remaining") | .ac // 0),
  ($v[$dt].populations[]? | select(.id=="remaining") | .an // 0),
  ($v[$dt].populations[]? | select(.id=="remaining") | .hemizygote_count // 0),
  ($v[$dt].populations[]? | select(.id=="remaining") | .homozygote_count // 0)
]
| @tsv
' "$input_file" >> "$output_file"

echo "✅ Output written to: $output_file"

