# master script to run the ortholog relationship data for all species
#
APPDIR=/home/rgddata/pipelines/OrthologRelationLoading
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

#load human-to-species orthologs (only direct orthologs)
#              "mouse" "rat" "dog" "bonobo" "squirrel" "chinchilla" "pig" "vervet" "molerat" )
SPECIES_LIST=( "2"     "3"   "6"   "5"      "7"        "4"          "9"   "13"     "14" )

# OLD: unportable: species common name can change
# SPECIES_LIST=( "mouse" "rat" "dog" "bonobo" "squirrel" "chinchilla" "pig" "vervet" "molerat" )

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
