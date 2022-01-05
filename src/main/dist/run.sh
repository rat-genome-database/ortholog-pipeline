#!/usr/bin/env bash
#
# ortholog loading pipeline wrapper script calling gradle-generated script
#
. /etc/profile
APPNAME=OrthologRelationLoading

APPDIR=/home/rgddata/pipelines/$APPNAME
cd $APPDIR

java -Dspring.config=$APPDIR/../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$APPDIR/properties/log4j2.xml \
    -jar lib/${APPNAME}.jar "$@" 2>&1 | tee run.log