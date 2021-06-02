#!/bin/bash
FDPS_SVN_USER=namekata
FDPS_SVN_HOST=v1.jmlab.jp
FDPS_SVN_WORKCOPY=fdps-svn
WORKDIR=work_for_upload
mkdir -p ${WORKDIR}; cd ${WORKDIR}

# [1] Checkout the latest SVN repository
svn checkout svn+ssh://${FDPS_SVN_USER}@${FDPS_SVN_HOST}//home/fdps/svn/fdps ${FDPS_SVN_WORKCOPY}
# Check this repository
# [1-1] Firstly build phantom-GRAPE library
#cd ${FDPS_SVN_WORKCOPY}/src/phantom_grape_x86/G5/newton/libpg5
#make
#cd -
## [1-2] Then, perform tests in sample/c++/nbody
#cd ${FDPS_SVN_WORKCOPY}/sample/c++/nbody
#./test.py
#if [ "$?" -ne 0 ]
#then
#   echo 'Failed to pass tests in sample/c++/nbody!'
#   echo 'Please resolve errors!'
#   exit 1
#fi
#cd -

# [2] Get a git clone of FDPS
git clone https://github.com/FDPS/FDPS.git 

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
git rm --force -r LICENSE
git rm --force -r README.md
rm -rf doc sample scripts src tests LICENSE README.md

# [3-3] Copy from a FDPS svn work copy
mkdir doc
cp -R ../${FDPS_SVN_WORKCOPY}/doc/doc_specs.pdf        doc/
cp -R ../${FDPS_SVN_WORKCOPY}/doc/doc_specs_en.pdf     doc/
cp ../${FDPS_SVN_WORKCOPY}/doc/doc_specs_cpp_ja.pdf    doc/
cp ../${FDPS_SVN_WORKCOPY}/doc/doc_specs_cpp_en.pdf    doc/
cp ../${FDPS_SVN_WORKCOPY}/doc/doc_specs_ftn_ja.pdf    doc/
cp ../${FDPS_SVN_WORKCOPY}/doc/doc_specs_ftn_en.pdf    doc/
cp -R ../${FDPS_SVN_WORKCOPY}/doc/doc_tutorial.pdf     doc/
cp -R ../${FDPS_SVN_WORKCOPY}/doc/doc_tutorial_e.pdf   doc/
cp ../${FDPS_SVN_WORKCOPY}/doc/doc_tutorial_cpp_ja.pdf doc/
cp ../${FDPS_SVN_WORKCOPY}/doc/doc_tutorial_cpp_en.pdf doc/
cp ../${FDPS_SVN_WORKCOPY}/doc/doc_tutorial_ftn_ja.pdf doc/
cp ../${FDPS_SVN_WORKCOPY}/doc/doc_tutorial_ftn_en.pdf doc/
cp ../${FDPS_SVN_WORKCOPY}/doc/doc_tutorial_c_ja.pdf   doc/
cp ../${FDPS_SVN_WORKCOPY}/doc/doc_tutorial_c_en.pdf   doc/
cp ../${FDPS_SVN_WORKCOPY}/doc/doc_SC15_handout.pdf    doc/
cp ../${FDPS_SVN_WORKCOPY}/doc/doc_SC15_handout.pptx   doc/
cp ../${FDPS_SVN_WORKCOPY}/LICENSE   ./
cp ../${FDPS_SVN_WORKCOPY}/README.md ./
svn export ../${FDPS_SVN_WORKCOPY}/src/      src
svn export ../${FDPS_SVN_WORKCOPY}/scripts/  scripts
svn export ../${FDPS_SVN_WORKCOPY}/sample/   sample
svn export ../${FDPS_SVN_WORKCOPY}/tests/    tests
# Delete unnecessary files
rm -rf ./src/fortran_interface-v1
rm -rf ./src/particle_mesh_multipole
rm -rf ./sample/c++/pmmm
rm -rf ./scripts/gen_ftn_if-v1.py

# [3-4] Label `add` to newly-created files
git add doc
git add src
git add scripts
git add sample
git add tests
git add LICENSE
git add README.md

# [4] Check status
git status
echo "If no problems, please do the followings:" 
echo "   $ git commit -m \"some message\"" 
echo "   $ git push" 
