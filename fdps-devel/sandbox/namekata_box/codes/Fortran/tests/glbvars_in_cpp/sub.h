#pragma once

//==============================
//* Class Decl.: Test object
//==============================
class Test_Object {
   private:
      // Variables
      int argc;
      char **argv;
      // Methods
   public:
      // Variables
      double var;
      // Methods
      void Initialize(int argc, char *argv[]);
};

extern double glbvar;
extern Test_Object glbobj;
