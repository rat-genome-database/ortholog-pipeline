#!/usr/bin/env bash
#
. /etc/profile
APPNAME=ortholog-pipeline

# exit with java's status, not tee's: without pipefail a failed run would look successful to the callers
# (f.e. runAll.sh aborts the per-species load only when the AGR load exits non-zero)
set -o pipefail

APPDIR=/home/rgddata/pipelines/$APPNAME
cd $APPDIR

java -Dspring.config=$APPDIR/../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$APPDIR/properties/log4j2.xml \
    -jar lib/${APPNAME}.jar "$@" 2>&1 | tee run.log