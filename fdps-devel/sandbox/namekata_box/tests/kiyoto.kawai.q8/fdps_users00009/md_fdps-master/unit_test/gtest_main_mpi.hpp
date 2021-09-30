//=======================================================================================
//  This is main function for unit test with google C++ test framework.
//=======================================================================================

//==========================================
// MPI initializer for google test
//  ref: http://www.parresianz.com/mpi/c++/mpi-unit-testing-googletests-cmake/
//==========================================
class MPIEnvironment :
    public ::testing::Environment
{
    public:
        virtual void SetUp(){
            char** argv;
            int argc = 0;
            PS::Initialize(argc, argv);

            //--- display total threads for FDPS
            if(PS::Comm::getRank() == 0){
                std::cerr << "Number of processes          : " << PS::Comm::getNumberOfProc()   << "\n"
                          << "Number of threads per process: " << PS::Comm::getNumberOfThread() << std::endl;
            }
        }

        virtual void TearDown(){
            PS::Finalize();
        }

        virtual ~MPIEnvironment() = default;
};

//==========================================
// main routine with MPI
//==========================================
int main(int argc, char* argv[]){
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::AddGlobalTestEnvironment(new MPIEnvironment);
    return RUN_ALL_TESTS();
}
