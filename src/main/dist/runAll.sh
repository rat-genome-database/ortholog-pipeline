# master script to run the ortholog relationship data for all species that are public in RGD
#  (SPECIES_TYPES.IS_SEARCHABLE=1)
# load human-to-species and species-to-human orthologs (only direct orthologs)
#
# Alliance orthologs are refreshed first because the per-species HCOP/NCBI cascade
# consults AGR_ORTHOLOGS as a priority-2 candidate (manual > Alliance > HCOP > NCBI).
# A stale or empty AGR_ORTHOLOGS table will cause the per-species run to abort.
#
APPDIR=/home/rgddata/pipelines/ortholog-pipeline
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

${APPDIR}/run.sh --agrOrthologs > ${APPDIR}/agrOrthologs.log || { echo "AGR ortholog load failed -- aborting"; exit 1; }

export SKIP_EMAIL_SUMMARY=1
$APPDIR/loadSpecies.sh all "$@"

#send summary email
EMAIL=mtutaj@mcw.edu
if [ "$SERVER" = "REED" ]; then
  EMAIL_LIST="mtutaj@mcw.edu slaulederkind@mcw.edu"
fi
mailx -s "[$SERVER] RGD orthologs OK" $EMAIL < ${APPDIR}/logs/status.log
