//=======================================================================================
//  This is main function for unit test with google C++ test framework.
//=======================================================================================

//==========================================
// main routine without MPI
//==========================================
int main(int argc, char* argv[]){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
