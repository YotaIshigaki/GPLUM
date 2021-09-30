== Updates since Code-4.3
* Bcast in initializePosDomain() is replaced by Allgather.
* New codes in ring/long_box/CodeIntegration4/ are merged.
* GetPosDomain3() added.
  This function is based on fdps/sample/c++/nbody-with-center/ring_calc.c

== Updates since CodeIntegration4.2
* GetPosDomain2(), which is a faster version of GetPosDomain(), is implemented.
  And the number of sampling points is increased.
