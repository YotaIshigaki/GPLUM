#!/bin/csh -f
#
# register spam on v1
#
# run this script at ~fdps/Mail/inbox as user fdps on v1.jmlab.jp
echo $1
cat $1 | nkf -m -w | bogofilter -s
cat $1 | nkf -m -w | bogofilter -s
cat $1 | nkf -m -w | bogofilter -v

