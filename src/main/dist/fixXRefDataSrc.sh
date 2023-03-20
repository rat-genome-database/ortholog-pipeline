# Diagnostics:
# 1) fix duplicates in RGD_ASSOCIATIONS.ASSOC_TYPE fields for associations of type 'weak_ortholog'
# 2) fix duplicates in GENETOGENE_RGD_ID_RLT.XREF_DATA_SET
# Reason: HCOP is emitting duplicate fields in its FTP file
# These duplicates for some genes cause erroneous ortholog assignments.
#
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
HOMEDIR=/home/rgddata/pipelines/ortholog-pipeline
EMAIL=mtutaj@mcw.edu

${HOMEDIR}/run.sh --fixXRefDataSet > fix.log
mailx -s "[$SERVER] RGD ortholog pipeline: source fix complete" $EMAIL < fix.log
