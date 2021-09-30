/*
 *  PZC Perfomance object
 *
 *	Copyright (c) 2014 PEZY Computing, K.K.
 */

#pragma once

#if defined(__PZC_KERNEL__)
	#include <pzc_builtin.h>
#endif

class	PZCPerformance
{
public:
	unsigned int	perf;		//!< freerun counter
	unsigned int	stall;		//!< stall counter
	unsigned int	wait;		//!< wait counter
	
	unsigned int    reserved[5];
	
public:
	PZCPerformance();
	~PZCPerformance();

	void Update( void ) volatile;
	void Clear( void );

	void Average( int num );

	PZCPerformance& operator  = ( const PZCPerformance& v );
	PZCPerformance& operator += ( const PZCPerformance& v );
	PZCPerformance& operator -= ( const PZCPerformance& v );

	unsigned int Perf ( void ) const;
	unsigned int Stall( void ) const;
	unsigned int Wait ( void ) const;

};

inline
unsigned int
PZCPerformance::Perf( void ) const
{
	return perf;
}

inline
unsigned int
PZCPerformance::Wait( void ) const
{
	return wait;
}

inline
unsigned int
PZCPerformance::Stall( void ) const
{
	return stall;
}

inline
void PZCPerformance::Clear( void )
{
	perf    = 0;
	stall   = 0;
	wait    = 0;
	reserved[0] = 0;
	reserved[1] = 0;
	reserved[2] = 0;
	reserved[3] = 0;
	reserved[4] = 0;
}


inline PZCPerformance::PZCPerformance( void )
{
	Clear();
}

inline PZCPerformance::~PZCPerformance( void )
{
}

inline
void
PZCPerformance::Update( void ) volatile
{
#if defined(__PZC_KERNEL__)
	register int _perf  = __builtin_pz_read_gpr(0x18);
	register int _stall = __builtin_pz_read_gpr(0x1a);
	register int _wait  = __builtin_pz_read_gpr(0x1e);

	perf    = _perf;
	stall   = _stall;
	wait    = _wait;
#endif
}

inline
PZCPerformance& PZCPerformance::operator  = ( const PZCPerformance& v )
{
	this->perf    = v.perf;
	this->stall   = v.stall;
	this->wait    = v.wait;

	return *this;
}

inline
void PZCPerformance::Average( int num )
{
	this->perf  /= num;
	this->stall /= num;
	this->wait  /= num;
}

inline
PZCPerformance& PZCPerformance::operator  += ( const PZCPerformance& v )
{
	this->perf    += v.perf;
	this->stall   += v.stall;
	this->wait    += v.wait;
	return *this;
}

inline
PZCPerformance& PZCPerformance::operator  -= ( const PZCPerformance& v )
{
	this->perf    -= v.perf;
	this->stall   -= v.stall;
	this->wait    -= v.wait;

	return *this;
}

inline PZCPerformance operator + ( const PZCPerformance& a, const PZCPerformance& b )
{
	PZCPerformance v = a;

	v += b;

	return v;
}

inline PZCPerformance operator - ( const PZCPerformance& a, const PZCPerformance& b )
{
	PZCPerformance v = a;

	v -= b;

	return v;
}


