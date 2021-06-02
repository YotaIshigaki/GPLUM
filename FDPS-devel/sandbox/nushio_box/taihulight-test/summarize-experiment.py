#!/usr/bin/env python3

import glob
import os
import re
import subprocess

def src_sortkey(fn):
    if fn == "Makefile":
        return (0,fn)
    elif re.search("\.h$", fn):
        return (1,fn)
    elif re.search("\.sh$", fn):
        return (9,fn)
    elif re.search("\.c$", fn) or  re.search("\.cpp$", fn):
        return (2,fn)
    else:
        return (-1,fn)

def is_source_file(fn):
    return src_sortkey(fn)[0] >= 0

def experiment_to_latex(log_fn):
    rootdir = re.sub("/log.txt", "/", log_fn)
    srcdir = rootdir + "src/"
    src_files = glob.glob(srcdir + "*")


    intro_tex = ""
    discuss_tex = ""
    with open(rootdir + "remark.tex", "r") as fp:
        src = fp.read().split("\n")
        title_tex = src[0]
        flag_discuss = False
        for l in src[1:]:
            if flag_discuss:
                discuss_tex += l + "\n"
            elif "%%%%"  in l:
                flag_discuss = True
            else:
                intro_tex += l + "\n"

    method_tex = ""
    for fn in sorted(src_files, key= src_sortkey):
        if not is_source_file(fn):
            continue

        with open(fn, "r") as fp:
            src_str = fp.read()
            method_tex += """
\\verb`{}`
\\begin{{code}}
{}
\\end{{code}}
""".format(fn.split("/")[-1],src_str)

    with open(log_fn, "r") as fp:
        src_str = fp.read()
        result_tex = """
Got the following results:

\\begin{{code}}
{}
\\end{{code}}
""".format(src_str)

    return "\n".join(
        ["\\subsection{" + title_tex + "}",
         "\\subsubsection{Introduction}",
         intro_tex,
         "\\subsubsection{Source code}\n We have compiled and executed the following source codes:",
         method_tex,
         "\\subsubsection{Results}",
         result_tex,
         "\\subsubsection{Discussion}",
         discuss_tex
     ])


for suffix in ["exercise", "practical"]:

    so,_ = subprocess.Popen("find src-{}/ -name log.txt".format(suffix), \
                            shell=True, stdout=subprocess.PIPE).communicate()

    tex = ""
    for fn in sorted(so.decode().split("\n")):
        if len(fn.strip()) <= 0:
            continue
        tex += experiment_to_latex(fn) + "\n"

    with open("content-{}.tex".format(suffix),"w") as fp:
        fp.write(tex)
