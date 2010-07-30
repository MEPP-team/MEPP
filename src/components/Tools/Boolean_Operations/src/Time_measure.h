#ifndef TIME_MEASURE_H
#define TIME_MEASURE_H

/*!
 * \file Time_measure.h
 * \brief An Object to measure the time
 * \author Cyril Leconte
 */

#ifdef _MSC_VER

#include <windows.h>

/*! \class Time_measure
 * \brief A timer
 *
 * This class creates a timer
 */
class Time_measure
{
private:
	/*! \brief clock frequency*/
	__int64 freq;
	/*! \brief initial time*/
	__int64 t0;

public:
	/*!
     *  \brief Constructor
	 *  
	 *  Constructor for Time_measure
     */
	Time_measure()
	{
		QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
	}

    /*!
     *  \brief Start the timer
	 *
	 *  This method resets the timer 
     */
	void Start()
	{
		QueryPerformanceCounter((LARGE_INTEGER*)&t0);
	}

    /*!
     *  \brief Get the elapsed time
	 *
	 *  Measures the elapsed time 
	 *
     *  \return the elapsed time in second
     */
	double GetDiff()
	{
		__int64 t1;
		QueryPerformanceCounter((LARGE_INTEGER*)&t1);
		return (double)(t1 - t0) / freq; 
	}
};
#else

#include <ctime>

/*! \class Time_measure
 * \brief A timer
 *
 * This class creates a timer
 */
class Time_measure
{
private:
	/*! \brief initial time*/
	double t0;

public:
	/*!
     *  \brief Constructor
	 *  
	 *  Constructor for Time_measure
     */
	Time_measure() {}

    /*!
     *  \brief Start the timer
	 *
	 *  This method resets the timer 
     */
	void Start() {t0 = clock();}

    /*!
     *  \brief Get the elapsed time
	 *
	 *  Measures the elapsed time 
	 *
     *  \return the elapsed time in second
     */
	double GetDiff()
	{
		return (double)((clock()-t0)/CLOCKS_PER_SEC);
	}
};

#endif // _MSC_VER

#endif // TIME_MEASURE_H
