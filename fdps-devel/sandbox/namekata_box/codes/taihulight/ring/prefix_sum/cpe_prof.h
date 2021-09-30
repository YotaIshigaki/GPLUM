#pragma once

// MPE cycle counts
static inline unsigned long rpcc_() {
  unsigned long rpcc;
  asm("rtc %0": "=r" (rpcc) : );
  return rpcc;
}

// CPE cycle counts
static inline unsigned long rtc_() {
  unsigned long rpcc;
  asm volatile("rcsr %0, 4":"=r"(rpcc));
  return rpcc;
}

// get CPE row id
void get_row_id_(int *row_id){
  int row;
  asm volatile("rcsr %0, 1" : "=r"(row));
  *row_id=row;
}

// get CPE column id
void get_col_id_(int *col_id){
  int col;
  asm volatile("rcsr %0, 2" : "=r"(col));
  *col_id=col;
}

// all CPE sync
#if 0
void sync_array_(){
  unsigned long sync_tmp;
  asm volatile(                  \
               "ldi %0, 0xff\n" \
               "sync %0\n" \
               "synr %0\n" \
               : \
               :"r"(sync_tmp):"memory");
}
#else
void sync_array_(){
  unsigned long sync_tmp;
  asm volatile(                  
	   "ldi %0, 0xff\n" 
	   "sync %0\n" 
	   "synr %0\n" 
	   :"=r"(sync_tmp) 
	   :
	   :"memory");
}
#endif

#define REG_PUTR(var, dst) asm volatile ("putr %0,%1\n"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1\n"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile ("getr %0\n":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0\n":"=r"(var))
#define NOP() asm volatile ("nop\n")
