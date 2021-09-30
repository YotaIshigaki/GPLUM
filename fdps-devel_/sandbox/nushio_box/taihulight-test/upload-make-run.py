#!/usr/bin/env python3

import datetime
import os
import re
import subprocess
import sys

remote_rootdir = "/home/export/base/nsccwuxi_riken/riken/online1/sandbox/nushio_box/sunway-test/"

local_workdir = os.getcwd()
local_srcdir = local_workdir + "/src/"
m = re.search("sunway-test\/(.*)", local_workdir)
uniquedir = m.group(1)
remote_srcdir = remote_rootdir + uniquedir + "/src/"

log_filename = "log.txt"

with open(log_filename, "w") as fp:
    fp.write("# {}\n".format(datetime.datetime.now()))

def shell(cmd):
    cmdmsg = "$ {}\n".format(cmd)
    sys.stderr.write(cmdmsg)
    with open(log_filename, "a") as fp:
        fp.write(cmdmsg)

    subprocess.call(cmd + " | tee -a " + log_filename, shell = True)

# generate the runner script

runscript_filename = local_srcdir + "run.sh"
with open(runscript_filename, "w") as fp:
    fp.write("""
cd {srcdir}
make && make run
    """.format(srcdir = remote_srcdir))
shell("chmod 755 " + runscript_filename)

shell("ssh sunway 'mkdir -p {}'".format(remote_srcdir))
shell("rsync -avz {} sunway:{}".format(local_srcdir, remote_srcdir))
shell("ssh sunway {}/run.sh".format(remote_srcdir))
