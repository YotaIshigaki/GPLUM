#!/bin/bash
  echo "Running on GDB on node `hostname`"
xterm -e gdb -x gdb_cmd --args $*
exit 0

