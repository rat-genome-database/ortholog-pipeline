# master script to run the ortholog relationship data for all species
#
APPDIR=/home/rgddata/pipelines/OrthologRelationLoading
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

#load human-to-species orthologs (only direct orthologs)
#  species name must be exactly as it is in SPECIES_TYPE.SHORT_NAME
SPECIES_LIST=( "mouse" "rat" "dog" "bonobo" "squirrel" "chinchilla" "pig" "vervet" "molerat" )

export SKIP_EMAIL_SUMMARY=1
for SPECIES in "${SPECIES_LIST[@]}"; do
    $APPDIR/loadSpecies.sh "$SPECIES" "$@"
done

#send summary email
EMAIL=mtutaj@mcw.edu
if [ "$SERVER" = "REED" ]; then
  EMAIL_LIST=mtutaj@mcw.edu,slaulederkind@mcw.edu
fi
mailx -s "[$SERVER] RGD orthologs OK" $EMAIL < ${APPDIR}/logs/status.log
