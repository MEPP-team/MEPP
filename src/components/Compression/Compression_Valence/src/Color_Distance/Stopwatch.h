/***************************************************************************
                                 Stopwatch.h
                             -------------------
    update               : 2003-04-25
    copyright            : (C) 2002-2003 by Michaël Roy
    email                : michaelroy@users.sourceforge.net
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef _STOPWATCH_
#define _STOPWATCH_

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include <iostream>
#include <iomanip>
#include <ctime>

//--
//
// Clock
//
//--
// Class Clock use in Stopwatch class (more convenient for output purpose).
class Clock
{

	//
	// Member variable
	//
	protected :

		clock_t time;
		
	//
	// Member functions
	//
	public :

		//
		// Constructors / Destructor
		//
		inline Clock() : time(0) {}
		inline Clock(const Clock& c) : time(c.time) {}
		inline Clock(const clock_t& t) : time(t) {}
		inline ~Clock() {}

		//
		// Time management
		//
		inline const clock_t& GetCurrentClock() { return time = clock(); }
		
		//
		// Accessors
		//
		inline clock_t& Time() { return time; }
		inline const clock_t& Time() const { return time; }

		//
		// Operators
		//
		inline operator clock_t&() { return time; }
		inline operator const clock_t&() const { return time; }

		inline Clock& operator=(const Clock& c) { time = c.time; return *this; }
		
		inline bool operator==(const Clock& c) const { return time == c.time; }
		inline bool operator!=(const Clock& c) const { return time != c.time; }
	
		inline Clock& operator+=(const Clock& c) { time += c.time; return *this; }
		inline Clock& operator-=(const Clock& c) { time -= c.time; return *this; }
		inline Clock& operator*=(const Clock& c) { time *= c.time; return *this; }
		inline Clock& operator/=(const Clock& c) { time /= c.time; return *this; }

		inline Clock operator+(const Clock& c) const { return Clock(*this) += c; }
		inline Clock operator-(const Clock& c) const { return Clock(*this) -= c; }
		inline Clock operator*(const Clock& c) const { return Clock(*this) *= c; }
		inline Clock operator/(const Clock& c) const { return Clock(*this) /= c; }
};

//--
//
// Output stream operator for the Clock class
//
//--
inline std::ostream& operator<<(std::ostream& os, const Clock& clock)
{
	return os<<clock.Time()/CLOCKS_PER_SEC<<'.'<<std::setw(2)<<std::setfill('0')<<(int)((float)clock.Time()/CLOCKS_PER_SEC*100.0)%100<<std::setfill(' ')<<"s";
}

//--
//
// Stopwatch
//
//--
// Implement a stopwatch in order to measure process time.
class Stopwatch
{
	
	//
	// Member variables
	//
	protected :
	
		bool on;     // Activity
		Clock now;   // Current time
		Clock last;  // Wall clock of start watch
		Clock inter; // Intermediate elapsed time
		Clock total; // Total elapsed time

	//
	// Member functions
	//
	public :
		
		// Default constructor
		inline Stopwatch() : on(false) {
		}
		
		// Destructor
		inline ~Stopwatch() {
		}
     
     		// Return activity of stopwatch
		inline bool IsRunning() const {
			return on;
		}
	
		// Stop stopwatch & initialize data
		inline Stopwatch& Reset() {
			on = false;
			inter.Time() = 0;
			total.Time() = 0;
			return *this;
		}

		// Active stopwatch
		inline Stopwatch& Start() { 
			if ( on == false )
			{
				last.GetCurrentClock();
				on = true; 
			}
			return *this;
		}
		
		// Unactive stopwatch
		inline Stopwatch& Stop() {
			if ( on == true )
			{
				now.GetCurrentClock();
				inter  = now;
				inter -= last;
				last   = now;
				total += inter;
				on     = false;
			}
			return *this;
		}


		// Return total time
		inline const Clock& Total() {
			if ( on == true )
			{	
				now.GetCurrentClock();
				inter  = now;
				inter -= last;
				last   = now;
				total += inter;
			}
			return total;
		}
	
		// Return intermediate time
		inline const Clock& Intermediate() {
			if ( on == true )
			{	
				now.GetCurrentClock();
				inter  = now;
				inter -= last;
				last   = now;
				total += inter;
			}
			return inter;
		}
};

#endif

#endif // _STOPWATCH_

