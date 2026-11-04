#!/bin/bash

# Example usage: ./shared/myhelpers_jauger/retrieve_from_research_pacs.sh 71451428 20260921
# Run only the local NRRD export after sorting has already completed:
# ./shared/myhelpers_jauger/retrieve_from_research_pacs.sh 71451428 20260921 STEP4_ONLY=1
# Backward-compatible environment-variable form:
# RESEARCH_PACS_STEP4_ONLY=1 ./shared/myhelpers_jauger/retrieve_from_research_pacs.sh 71451428 20260921

usage() {
  echo "Usage: $0 <MRN> <StudyDate> [STEP4_ONLY=1]" >&2
  echo "  StudyDate must use YYYYMMDD format." >&2
  echo "  STEP4_ONLY=1 skips PACS retrieval and Steps 1-3." >&2
  echo "  RESEARCH_PACS_STEP4_ONLY=1 remains supported for backward compatibility." >&2
}

if [[ $# -lt 2 || $# -gt 3 ]]; then
  usage
  exit 2
fi

MRN="$1"
StudyDate="$2"
STEP4_ONLY="${RESEARCH_PACS_STEP4_ONLY:-0}"

if [[ $# -eq 3 ]]; then
  if [[ $3 != "STEP4_ONLY=1" ]]; then
    echo "Error: unrecognized third argument '$3'; expected STEP4_ONLY=1." >&2
    usage
    exit 2
  fi
  STEP4_ONLY=1
fi

if [[ ! ${StudyDate} =~ ^[0-9]{8}$ ]]; then
  echo "Error: StudyDate must be an 8-digit date in YYYYMMDD format." >&2
  usage
  exit 2
fi

if [[ ${STEP4_ONLY} != "0" && ${STEP4_ONLY} != "1" ]]; then
  echo "Error: RESEARCH_PACS_STEP4_ONLY must be 0 or 1." >&2
  exit 2
fi

HOST_DIRECTORY="$(pwd)"
HOST_UID="$(id -u)"
HOST_GID="$(id -g)"

if [[ ${STEP4_ONLY} == "0" ]]; then
  # Step 0: Pull data from Research Synapse.
  sudo docker run --rm -it -u "${HOST_UID}:${HOST_GID}" \
  --volume "${HOST_DIRECTORY}:/data" crl/dicom-tools:latest retrieve_dicoms.py \
  --outputDir "/data/${MRN}" \
  --aec SYNAPSERESEARCH --aet PACSDCM --namednode 10.20.2.28 --modality MR \
  --subjectID "${MRN}" --studyDate "${StudyDate}"

  # Step 1: Uncompress the 4D enhanced DICOMs to separate 3D dicoms.
  sudo docker run --rm -u "${HOST_UID}:${HOST_GID}" \
    -v "${HOST_DIRECTORY}:/data" crl/dicom-tools \
    uncompress_dicoms.py "./${MRN}" "./${MRN}-1-uncompressed"

  # Step 2: Sort the uncompressed 3D DICOMs into each sequence.
  sudo docker run --rm -u "${HOST_UID}:${HOST_GID}" \
    -v "${HOST_DIRECTORY}:/data" crl/dicom-tools \
    sort_dicoms.py "./${MRN}-1-uncompressed" "./${MRN}-2-sorted"

  # Step 3: Convert each 3D DICOM sequence into a single 4D NIfTI file
  # with its JSON metadata file.
  sudo docker run --rm -u "${HOST_UID}:${HOST_GID}" \
    -v "${HOST_DIRECTORY}:/data" crl/dicom-tools \
    dicom_tree_to_nifti.py "./${MRN}-2-sorted" "${MRN}-3-converted"
else
  echo "STEP4_ONLY=1: skipping Research Synapse retrieval and Steps 1-3."
fi

validate_nrrd_file() {
  local nrrd_file="$1"
  local first_line
  local data_spec
  local data_path
  local listed_file
  local listed_count=0
  local listed_missing=0

  if ! IFS= read -r first_line < "${nrrd_file}" || [[ ${first_line} != NRRD* ]]; then
    echo "  Invalid NRRD header: ${nrrd_file}" >&2
    return 1
  fi

  data_spec="$(awk '
    NR > 1 && $0 ~ /^\r?$/ { exit }
    tolower($0) ~ /^[[:space:]]*data[[:space:]]+file[[:space:]]*:/ {
      line=$0
      sub(/^[^:]*:[[:space:]]*/, "", line)
      sub(/\r$/, "", line)
      print line
      exit
    }
  ' "${nrrd_file}")"

  if [[ -z ${data_spec} ]]; then
    if [[ ${nrrd_file,,} == *.nhdr ]]; then
      echo "  Detached NRRD header has no data file field: ${nrrd_file}" >&2
      return 1
    fi
    return 0
  fi

  if [[ ${data_spec^^} == LIST* ]]; then
    while IFS= read -r listed_file; do
      [[ -z ${listed_file} ]] && continue
      listed_count=$((listed_count + 1))
      data_path="$(dirname "${nrrd_file}")/${listed_file}"
      if [[ ! -f ${data_path} ]]; then
        echo "  Missing detached NRRD data file: ${data_path}" >&2
        listed_missing=$((listed_missing + 1))
      fi
    done < <(awk '
      capture {
        line=$0
        sub(/\r$/, "", line)
        if (length(line) > 0) print line
      }
      tolower($0) ~ /^[[:space:]]*data[[:space:]]+file[[:space:]]*:[[:space:]]*list/ {
        capture=1
      }
    ' "${nrrd_file}")
    [[ ${listed_count} -gt 0 && ${listed_missing} -eq 0 ]]
    return
  fi

  if [[ ${data_spec} == *%* ]]; then
    echo "  Unsupported detached NRRD filename sequence: ${data_spec}" >&2
    return 1
  fi

  if [[ ${data_spec} = /* ]]; then
    data_path="${data_spec}"
  else
    data_path="$(dirname "${nrrd_file}")/${data_spec}"
  fi
  if [[ ! -f ${data_path} ]]; then
    echo "  Missing detached NRRD data file: ${data_path}" >&2
    return 1
  fi
}

validate_nrrd_directory() {
  local output_directory="$1"
  local nrrd_file
  local found=0
  local invalid=0

  while IFS= read -r -d '' nrrd_file; do
    found=$((found + 1))
    if ! validate_nrrd_file "${nrrd_file}"; then
      invalid=$((invalid + 1))
    fi
  done < <(
    find "${output_directory}" -mindepth 1 -maxdepth 1 -type f \
      \( -iname '*.nrrd' -o -iname '*.nhdr' \) -print0
  )

  if [[ ${found} -eq 0 ]]; then
    echo "  No .nrrd or .nhdr output was produced in ${output_directory}" >&2
    return 1
  fi
  [[ ${invalid} -eq 0 ]]
}

# Step 4: Export NRRD images directly from each sorted, uncompressed DICOM
# series. These files support comparison of PACS-derived geometry and pixel
# data with NRRDs saved natively by the motion-control-stack FIRE-server.
SORTED_DIRECTORY="${HOST_DIRECTORY}/${MRN}-2-sorted"
NRRD_DIRECTORY="${HOST_DIRECTORY}/${MRN}-4-nrrd"

if [[ ! -d ${SORTED_DIRECTORY} ]]; then
  echo "Error: sorted DICOM directory does not exist: ${SORTED_DIRECTORY}" >&2
  exit 1
fi
mkdir -p "${NRRD_DIRECTORY}"

successful_series=0
failed_series=0
discovered_series=0

while IFS= read -r -d '' series_directory; do
  has_child_directory=0
  has_file=0

  while IFS= read -r -d '' _child_directory; do
    has_child_directory=1
    break
  done < <(find "${series_directory}" -mindepth 1 -maxdepth 1 -type d -print0)
  [[ ${has_child_directory} -eq 1 ]] && continue

  while IFS= read -r -d '' _series_file; do
    has_file=1
    break
  done < <(find "${series_directory}" -mindepth 1 -maxdepth 1 -type f -print0)
  [[ ${has_file} -eq 0 ]] && continue

  discovered_series=$((discovered_series + 1))
  if [[ ${series_directory} == "${SORTED_DIRECTORY}" ]]; then
    relative_directory="."
    output_directory="${NRRD_DIRECTORY}"
    container_input="/data/${MRN}-2-sorted"
    container_output="/data/${MRN}-4-nrrd"
  else
    relative_directory="${series_directory#"${SORTED_DIRECTORY}/"}"
    output_directory="${NRRD_DIRECTORY}/${relative_directory}"
    container_input="/data/${MRN}-2-sorted/${relative_directory}"
    container_output="/data/${MRN}-4-nrrd/${relative_directory}"
  fi
  mkdir -p "${output_directory}"

  echo "Step 4: converting series: ${relative_directory}"
  echo "  Output directory: ${output_directory}"
  if sudo docker run --rm -u "${HOST_UID}:${HOST_GID}" \
      -v "${HOST_DIRECTORY}:/data" crl/dicom-tools:latest \
      dcm2niix -e y -z n -b n -o "${container_output}" "${container_input}"; then
    if validate_nrrd_directory "${output_directory}"; then
      echo "  SUCCESS: ${relative_directory}"
      successful_series=$((successful_series + 1))
    else
      echo "  FAILED validation: ${relative_directory}" >&2
      failed_series=$((failed_series + 1))
    fi
  else
    echo "  FAILED conversion: ${relative_directory}" >&2
    failed_series=$((failed_series + 1))
  fi
done < <(find "${SORTED_DIRECTORY}" -type d -print0 | sort -z)

if [[ ${discovered_series} -eq 0 ]]; then
  echo "Step 4: no leaf directories containing files were found." >&2
  exit 1
fi

echo "Step 4 complete: ${successful_series} successful, ${failed_series} failed."
if [[ ${failed_series} -ne 0 ]]; then
  exit 1
fi

exit 0

