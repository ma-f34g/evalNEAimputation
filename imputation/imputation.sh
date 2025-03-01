#!/bin/bash

BEAGLE_VERSION="17Dec24.224"
BASE_URL="http://faculty.washington.edu/browning/beagle"

BEAGLE_JAR="beagle.${BEAGLE_VERSION}.jar"
REF_DATA="../data/reference_panel.vcf.gz"
REF_PANEL="../data/reference_panel.bref3"

set -ex

for jar in "${BEAGLE_JAR}"; do
	if [ ! -f "${jar}" ]; then
		echo "Downloading ${jar}"
		wget "${BASE_URL}/${jar}"
	fi
done

bref3 $REF_DATA > $REF_PANEL
for TEST_DATA in ../data/damaged_test_samples_*.vcf.gz; do
	DAMAGE_RATE=$(basename "${TEST_DATA}" | grep -o '[0-9]\.[0-9]')
	REF_OUTPUT="../data/imputed_${DAMAGE_RATE}"

	echo "*** Processing file with damage rate ${DAMAGE_RATE} ***"
	bcftools view ${TEST_DATA} | sed 's/0|0/\.\/\./g' | bcftools view -Oz -o ${TEST_DATA}_cleaned.gz
	java -Xmx250g -jar "${BEAGLE_JAR}" ref="${REF_PANEL}" gt="${TEST_DATA}_cleaned.gz" out="${REF_OUTPUT}"
	tabix ${REF_OUTPUT}.vcf.gz
	bcftools filter -i 'INFO/IMP=1' ${REF_OUTPUT}.vcf.gz -Oz -o ${REF_OUTPUT}_filtered.vcf.gz
done
