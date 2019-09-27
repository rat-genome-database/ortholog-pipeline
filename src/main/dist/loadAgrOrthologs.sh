# Load orthologs from AGR
#
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
HOMEDIR=/home/rgddata/pipelines/OrthologRelationLoading
EMAIL=mtutaj@mcw.edu

${HOMEDIR}/run.sh --agrOrthologs > agrOrthologs.log

mailx -s "[$SERVER] RGD ortholog pipeline: AGR orthologs loaded" $EMAIL < agrOrthologs.log
