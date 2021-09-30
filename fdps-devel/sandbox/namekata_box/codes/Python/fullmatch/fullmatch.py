#!/usr/bin/env /usr/bin/python3.4
# -*- coding: utf-8 -*-
#=================================
#   Module import
#=================================
try:
    import os
    import sys
    import re
    import struct
    import math
except ImportError:
    print("Module os,sys,struct,math are not found.")
    quit()

#=================================
#   Function
#=================================
def __fullmatch(pattern, string, flags=0):
    """Emulate python-3.4 re.fullmatch(pattern, string, flags=0)."""
    m = re.match("(?:" + pattern + r")\Z", string, flags=flags)
    return m

#=================================
#   Main
#=================================
string  = "!$fdps FP,EPI,EPJ,Force"
pattern = r"^!\$fdps.*"
m = __fullmatch(pattern,string)
if (m):
    print(m)
else:
    print(m)
