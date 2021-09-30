== Updates since Code-4.4
* exchangeParticle4() added.
* Information output to fout_debug is increased.
* A simple load balancer (experimentally phase) is added.
* To identify the host causing "Unknown exception",
  a new macro GEN_MORTON_KEY_SEQ_EXEC_FOR_DEBUG and
  the corresponding code is added.

== Updates since Code-4.3
* Bcast in initializePosDomain() is replaced by Allgather.
* New codes in ring/long_box/CodeIntegration4/ are merged.
* GetPosDomain3() added.
  This function is based on fdps/sample/c++/nbody-with-center/ring_calc.c

== Updates since CodeIntegration4.2
* GetPosDomain2(), which is a faster version of GetPosDomain(), is implemented.
  And the number of sampling points is increased.
