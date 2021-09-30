#!/bin/bash
  echo "Running on GDB on node `hostname`"
xterm +bdc +cm -e gdb -x cmd_gdb --args $*
#xterm +bdc +cm -e ddd& -x cmd_gdb --args $*
exit 0
