extern "C" {

int call_it (int (*func)(int), int arg) {
       return func (arg);
}

} // END of extern "C"
