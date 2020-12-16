/* ./GRAPEEmulator.c */
void InitializeGRAPEEmulator(void);
void g5_open_emu(void);
void g5_close_emu(void);
int g5_get_number_of_pipelines_emu(void);
int g5_get_number_of_boards_emu(void);
void g5_set_range_emu(const double size1, const double size2, const double minimum_mass);
void g5_set_n_emu(const int number);
void g5_set_xmj_emu(const int offset, const int number, double Pos[restrict][3], double Mass[restrict]);
void g5_set_ip_emu(const int number_of_pipelines, double Pos[restrict][3], double Eps[restrict], double Hsml[restrict]);
void g5_run_emu(void);
void g5_get_force_emu(const int number_of_pipelines, double Acc[restrict][3], double Pot[restrict]);
