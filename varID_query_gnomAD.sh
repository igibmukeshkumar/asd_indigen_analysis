#!/bin/bash

# Function to display usage
usage() {
  echo "Usage: $0 --vl <variant_list.txt> --data-type <exome|genome|joint>"
  echo "  --vl        : Path to the variant list file (e.g., variants.txt)"
  echo "  --data-type : Choose data type (exome, genome, joint)"
  echo "  --help      : Display this help message"
  exit 1
}

# Parse command-line options
while [[ "$1" != "" ]]; do
  case "$1" in
    --vl)
      shift
      VARIANT_LIST_FILE="$1"
      ;;
    --data-type)
      shift
      DATA_TYPE="$1"
      ;;
    --help)
      usage
      ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
  esac
  shift
done

# Input validation
if [ -z "$VARIANT_LIST_FILE" ]; then
  echo "Error: Variant list file is required (--vl)"
  usage
fi

if [ -z "$DATA_TYPE" ] || [[ "$DATA_TYPE" != "exome" && "$DATA_TYPE" != "genome" && "$DATA_TYPE" != "joint" ]]; then
  echo "Error: Invalid data type. Please specify exome, genome, or joint (--data-type)"
  usage
fi

if [ ! -f "$VARIANT_LIST_FILE" ]; then
  echo "File not found: $VARIANT_LIST_FILE"
  exit 1
fi

# Make output directory
mkdir -p varStore

# Start processing
while read -r VARIANT; do
  [ -z "$VARIANT" ] && continue
  echo "🔍 Querying gnomAD for variant: $VARIANT"

  # Set fields dynamically based on data type
  if [ "$DATA_TYPE" == "joint" ]; then
    FIELDS=$(cat <<EOF
      $DATA_TYPE {
        ac
        an
        hemizygote_count
        homozygote_count
        populations {
          id
          ac
          an
          hemizygote_count
          homozygote_count
        }
      }
EOF
)
  else
    FIELDS=$(cat <<EOF
      $DATA_TYPE {
        ac
        an
        af
        ac_hemi
        ac_hom
        populations {
          id
          ac
          an
          ac_hemi
          ac_hom
        }
      }
EOF
)
  fi

  # Query gnomAD
  RESP=$(curl -s https://gnomad.broadinstitute.org/api \
    -H 'Content-Type: application/graphql; charset=utf-8' \
    --data-binary @- <<EOF
query VariantInfo {
  variant(variantId: "$VARIANT", dataset: gnomad_r4) {
    variant_id
    rsids
    chrom
    pos
    ref
    alt
    transcript_consequences {
      gene_symbol
      transcript_id
      transcript_version
      hgvsc
      hgvsp
      lof
      major_consequence
    }
    $FIELDS
  }
}
EOF
)

  # Save only if variant was found
  if echo "$RESP" | grep -q '"variant":null'; then
    echo "❌ Variant not found in gnomAD: $VARIANT"
  else
    echo "$RESP" > "varStore/${VARIANT//\//_}.json"
    echo "✅ Saved: varStore/${VARIANT//\//_}.json"
  fi

done < "$VARIANT_LIST_FILE"

echo "🎉 All found variants saved in varStore/"
