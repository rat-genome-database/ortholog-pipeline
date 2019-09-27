#!/usr/bin/env bash
#
# ortholog loading pipeline wrapper script calling gradle-generated script
#
. /etc/profile
APPNAME=OrthologRelationLoading

APPDIR=/home/rgddata/pipelines/$APPNAME
cd $APPDIR

java -Dspring.config=$APPDIR/../properties/default_db2.xml \
    -Dlog4j.configuration=file://$APPDIR/properties/log4j.properties \
    -jar lib/${APPNAME}.jar "$@" 2>&1 | tee run.log