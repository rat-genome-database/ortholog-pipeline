# Script to run the ortholog relationship data loading pipeline for given species
# human==>species and species==>human orthologs will be loaded
#
# example: ./loadSpecies.sh dog
#
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
HOMEDIR=/home/rgddata/pipelines/OrthologRelationLoading

${HOMEDIR}/run.sh --species $1 > ${HOMEDIR}/$1.log
