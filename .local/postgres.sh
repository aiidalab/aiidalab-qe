#!/bin/bash

PGBIN=/usr/lib/postgresql/10/bin

# -w waits until server is up
PSQL_START_CMD="${PGBIN}/pg_ctl --timeout=180 -w -D ${AIIDALAB_HOME}/.postgresql -l ${AIIDALAB_HOME}/.postgresql/logfile start"
PSQL_STOP_CMD="${PGBIN}/pg_ctl -w -D ${AIIDALAB_HOME}/.postgresql stop"
PSQL_STATUS_CMD="${PGBIN}/pg_ctl -D ${AIIDALAB_HOME}/.postgresql status"

# helper function to start psql and wait for it
#TIMEOUT=20
#until psql -h localhost template1 -c ";" || [ $TIMEOUT -eq 0 ]; do
#   echo ">>>>>>>>> Waiting for postgres server, $((TIMEOUT--)) remaining attempts..."
#   tail -n 50 ${AIIDALAB_HOME}/.postgresql/logfile
#   sleep 1
#done

function psql_start {
    # make DB directory, if not existent
    if [ ! -d ${AIIDALAB_HOME}/.postgresql ]; then
       mkdir ${AIIDALAB_HOME}/.postgresql
       ${PGBIN}/initdb -D ${AIIDALAB_HOME}/.postgresql
       echo "unix_socket_directories = '${AIIDALAB_HOME}/.postgresql'" >> ${AIIDALAB_HOME}/.postgresql/postgresql.conf
       ${PSQL_START_CMD}
       psql -h localhost -d template1 -c "CREATE USER aiida WITH PASSWORD 'aiida_db_passwd';"
       psql -h localhost -d template1 -c "CREATE DATABASE aiidadb OWNER aiida;"
       psql -h localhost -d template1 -c "GRANT ALL PRIVILEGES ON DATABASE aiidadb to aiida;"

    # else don't
    else

        # stores return value in $?
        running=true
        ${PSQL_STATUS_CMD} || running=false

        # Postgresql was probably not shutdown properly. Cleaning up the mess...
        if ! $running ; then
           echo "" > ${AIIDALAB_HOME}/.postgresql/logfile # empty log files
           rm -vf ${AIIDALAB_HOME}/.postgresql/.s.PGSQL.5432
           rm -vf ${AIIDALAB_HOME}/.postgresql/.s.PGSQL.5432.lock
           rm -vf ${AIIDALAB_HOME}/.postgresql/postmaster.pid
           ${PSQL_START_CMD}
       fi
    fi
}

function psql_stop {
    $PSQL_STOP_CMD
}
