# master script to run the ortholog relationship data for all species that are public in RGD
#  (SPECIES_TYPES.IS_SEARCHABLE=1)
# load human-to-species and species-to-human orthologs (only direct orthologs)
#
APPDIR=/home/rgddata/pipelines/ortholog-pipeline
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

export SKIP_EMAIL_SUMMARY=1
$APPDIR/loadSpecies.sh all "$@"

#send summary email
EMAIL=mtutaj@mcw.edu
if [ "$SERVER" = "REED" ]; then
  EMAIL_LIST=mtutaj@mcw.edu,slaulederkind@mcw.edu
fi
mailx -s "[$SERVER] RGD orthologs OK" $EMAIL < ${APPDIR}/logs/status.log
