#!/bin/bash

for op in faddd fsubd fmuld fmaxd fmind frsqrtad frcpad fmaddd fnmaddd fmsubd fnmsubd fmovd; do
  echo ${op}
  grep ${op} $1 | wc -l
done
