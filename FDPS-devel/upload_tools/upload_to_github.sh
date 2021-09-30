#!/bin/bash
FDPS_WORKCOPY=/home/masaki/project/git/fdps-devel/
WORKDIR=work_for_upload
mkdir -p ${WORKDIR}; cd ${WORKDIR}

# [2] Get a git clone of FDPS
#git clone https://github.com/FDPS/FDPS.git 
git clone git@github.com:FDPS/FDPS.git

# [3] Update local git clone using a FDPS svn work copy
cd FDPS

# [3-1] Set the git username and email for this repository
git config user.name "FDPS official" 
git config user.email "fdps_official@mail.jmlab.jp" 

# [3-2] Clear
git rm --force -r doc
git rm --force -r sample
git rm --force -r scripts
git rm --force -r src
git rm --force -r tests
git rm --force -r pikg
git rm --force -r LICENSE
git rm --force -r README.md
rm -rf doc sample scripts src tests pikg LICENSE README.md

# [3-3] Copy from a FDPS svn work copy
mkdir doc
cp -d ${FDPS_WORKCOPY}/doc/doc_specs.pdf           ./doc/
cp -d ${FDPS_WORKCOPY}/doc/doc_specs_cpp_ja.pdf    ./doc/
cp -d ${FDPS_WORKCOPY}/doc/doc_specs_en.pdf        ./doc/
cp -d ${FDPS_WORKCOPY}/doc/doc_specs_cpp_en.pdf    ./doc/
cp -d ${FDPS_WORKCOPY}/doc/doc_specs_ftn_ja.pdf    ./doc/
cp -d ${FDPS_WORKCOPY}/doc/doc_specs_ftn_en.pdf    ./doc/
cp -d ${FDPS_WORKCOPY}/doc/doc_tutorial.pdf        ./doc/
cp -d ${FDPS_WORKCOPY}/doc/doc_tutorial_e.pdf      ./doc/
cp -d ${FDPS_WORKCOPY}/doc/doc_tutorial_cpp_ja.pdf ./doc/
cp -d ${FDPS_WORKCOPY}/doc/doc_tutorial_cpp_en.pdf ./doc/
cp -d ${FDPS_WORKCOPY}/doc/doc_tutorial_ftn_ja.pdf ./doc/
cp -d ${FDPS_WORKCOPY}/doc/doc_tutorial_ftn_en.pdf ./doc/
cp -d ${FDPS_WORKCOPY}/doc/doc_tutorial_c_ja.pdf   ./doc/
cp -d ${FDPS_WORKCOPY}/doc/doc_tutorial_c_en.pdf   ./doc/
cp -d ${FDPS_WORKCOPY}/doc/doc_SC15_handout.pdf    ./doc/
cp -d ${FDPS_WORKCOPY}/doc/doc_SC15_handout.pptx   ./doc/

cp ${FDPS_WORKCOPY}/LICENSE   ./
cp ${FDPS_WORKCOPY}/README.md ./

cp -r ${FDPS_WORKCOPY}/src/ ./
rm -rf ./src/fortran_interface-v1
rm -rf ./src/particle_mesh_multipole

cp -r ${FDPS_WORKCOPY}/scripts/ ./
rm -rf ./scripts/gen_ftn_if-v1.py

cp -r ${FDPS_WORKCOPY}/sample/ ./
rm -rf ./sample/c++/pmmm

cp -r ${FDPS_WORKCOPY}/tests/ ./

cp -r ${FDPS_WORKCOPY}/pikg/ ./

# Delete unnecessary files
rm -rf ./src/fortran_interface-v1
rm -rf ./src/particle_mesh_multipole
rm -rf ./sample/c++/pmmm
rm -rf ./scripts/gen_ftn_if-v1.py
find ./sample/ -name "*.out" -type f | xargs rm

# [3-4] Label `add` to newly-created files
git add doc
git add src
git add scripts
git add sample
git add tests
git add pikg
git add LICENSE
git add README.md

# [4] Check status
git status
echo "If no problems, please do the followings:" 
echo "   $ git commit -m \"some message\"" 
echo "   $ git push" 

