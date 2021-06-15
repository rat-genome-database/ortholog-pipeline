# Script to run the ortholog relationship data loading pipeline for given species
# human==>species and species==>human orthologs will be loaded
#
# example: ./loadSpecies.sh dog
#
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
HOMEDIR=/home/rgddata/pipelines/OrthologRelationLoading
EMAIL=mtutaj@mcw.edu
if [ "$SERVER" = "REED" ]; then
  EMAIL_LIST=mtutaj@mcw.edu,slaulederkind@mcw.edu
fi

${HOMEDIR}/run.sh --species $1 > ${HOMEDIR}/$1.log
if [ $SKIP_EMAIL_SUMMARY -ne 1 ]; then
  mailx -s "[$SERVER] RGD orthologs OK for $1" $EMAIL < ${HOMEDIR}/$1.log
fi