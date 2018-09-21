#!/usr/bin/env bash
#
# ortholog loading pipeline wrapper script calling gradle-generated script
#
. /etc/profile
APPNAME=OrthologRelationLoading

APPDIR=/home/rgddata/pipelines/$APPNAME
cd $APPDIR
pwd
DB_OPTS="-Dspring.config=$APPDIR/../properties/default_db.xml"
LOG4J_OPTS="-Dlog4j.configuration=file://$APPDIR/properties/log4j.properties"
declare -x "ORTHOLOG_RELATION_LOADING_OPTS=$DB_OPTS $LOG4J_OPTS"
bin/$APPNAME "$@" 2>&1 | tee run.log