/****************************************************************
 * Histogram
 *
 * Objects of this class store general histograms. The class is
 * implemented as template, so the type can be choosen.
 *
 * IMPORTANT: the parameter 'discrete' of each object, given at
 * the constructor, tells the object if it stores discrete values
 * like a graylevel histogram, or continuous values, like the
 * real values of the range [0-1]. This makes a difference in the
 * stepsize:
 *
 * discrete:   stepsize = (max - min  +  1 ) / bincount
 *                  e.g.  (255 -  0   +  1 ) / 256       = 1
 * continuous: stepsize = (max - min ) / bincount
 *                  e.g.  (1.0 - 0.0 ) / 4               = 0.25
 *
 * A graylevel histogram which uses (wrongly) the continous form,
 * would result in a stepsize of 255/256, which is obviously
 * wrong.
 *
 * --------------------------------------------------------------
 *
 * Some operations make only sense on the cumulative histogram.
 * A histogram may be transformed into a cumulative histogram by
 * applying the toCumulative() method. However, the histogram
 * object does not store the fact that the call has been made,
 * so later it will note be aware of this fact. It is the user
 * who must manage this transitions.
 *
 * Christian Wolf
 * wolf@rfv.insa-lyon.fr
 ****************************************************************/

#ifndef _WOLF_HISTOGRAM_H_
#define _WOLF_HISTOGRAM_H_

#include <iostream>

// From Vincent Malleron
#include "../gnuplot/gnuplot_interface.h"

using namespace std;

template <class T>
class Histogram {

	public:

		// The datatype used in the histogram
		typedef T datatype;

		// Constructors
		Histogram (const Histogram<T> & other);
		Histogram (int s, T min, T max, bool isdiscret);
		//Histogram (const DataSerie &s);

		// Destructor
		~Histogram () { if (bins!=NULL) delete bins; }

		// Assignment operator
		Histogram<T> & operator = (const Histogram<T> & other);

		// Access to the bins
		int value2index(T val);
		T & operator [] (int index) { return bins[index]; }
		T & operator () (int index) { return bins[index]; }
		T & fromValue (T val) 		{ return bins[value2index(val)]; }

		// Add and remove a value - i.e. increase the histogram by one
		void add (T val);
		void remove (T val);

		// Add and remove a value - i.e. increase the histogram by a weight
		void add (T val, T weight);

		// - operator: L2 norm
		double operator - (const Histogram<T> & other);

		// standard operations
		void normalize();
		void clear();
		T sum () const;
		void stretch() ;

		// Statistical operations
		double getMean();
		double getVariance();

		// Create statistical distributions
		void createUniformDistributionFromValue (T val);
		void createUniformDistributionFromMass (T mass);

		// Transform the histogram into a cumulative histogram
		void toCumulative();
		void toCumulativeFromRight();

		void toLogScale();

        void truncFromRightStop(double percentage);

		// Operations on cumulative histograms
		float kolmogorovSmirnovTest(Histogram<T> &other) const;

		// Difference histogram
		Histogram<T> * diffHist (const Histogram<T> &other);

		// The earth mover's distance between two histograms
		double EMD (const Histogram<T> &other) const;

		inline int getThresholdValueFisher (int kmin, int kmax)
		{
			return getThresholdValueFisher (kmin, kmax, NULL, NULL);
		}
		int getThresholdValueFisher (int kmin, int kmax, float *out_peakleft, float *out_peakright);
		inline void getThresholdValueFisher2th (int kmin, int kmax, int &out_k1, int &out_k2);

		double getIntraclassVariance (int t);

		// Get the borders
		inline int firstNonZero();
		inline int lastNonZero();
		inline int getSize() { return _size; }

		// Get the index of the maximum
		int argmax();

		// Get the bin index of this value
		int getIndexOf(T val) { return (int) ((double) (val - min) / step) / 1; };

		// Get the value of the center of the given bin
		T getBinCenter(int bin);
		T getStepSize() 								{ return step; }
		T getMin() 										{ return min; }
		T getMax() 										{ return max; }

		int getIndexPercLeft (double perc);
		int getIndexPercRight (double perc);
		void cutOutliers (double leftperc, double rightperc);
		void removeFalsePeaks();

		// Filtering (1D convolution)
		//void filter (DataSerie &fm);

		// A special feature for histograms which store grayvalues:
		unsigned char getContrast ();

		// OUTPUT
		void matlabOut(ostream & os, const char *s);

        void gnuplotOut (const char *title, const char *Xlabel, const char *Ylabel, bool log_scale=false);

	private:

	public:
		int _size;
		T min, max;
		double step;
		bool discret;

		T *bins;
};

// Global functions
template <class T> ostream & operator << (ostream & os, Histogram<T> h);

// *************************************************************
// Constructor - empty Histogram
// *************************************************************

template <class T>
inline Histogram<T>::Histogram (int s, T mi, T ma, bool d) {
   	_size = s;
   	min = (T) mi;
	max = (T) ma;
	discret = d;

	// The discret case
	if (discret)
		step = ((double)((double)max-(double)min+1.0))/ (double) _size;
	else
		step = ((double)((double)max-(double)min))/ (double) _size;

   	bins = new T[s];
   	clear();
}


// *************************************************************
// Constructor - Convert from a DataSerie
// *************************************************************
/*
template <class T>
inline Histogram<T>::Histogram (const DataSerie &s) {

    _size = s.size();
    min = 0;
    max = 255;
    discret = true;
    step=1;

    bins = new T [_size];
    clear();

    for (int x=0; x<_size; ++x) {
    	bins[x] = (T) s[x];
    	if (min>bins[x])
    		min=bins[x];
		if (max<bins[x])
    		max=bins[x];
    }
}
*/

// *************************************************************
// Copy constructor
// *************************************************************

template <class T>
inline Histogram<T>::Histogram (const Histogram<T> & other) {

	_size = other._size;
	min = other.min;
	max = other.max;
	step = other.step;
	discret = other.discret;
	bins = new T [_size];

	for(int i=0; i<_size; ++i)
		bins[i] = other.bins[i];
}


// *************************************************************
// Assignment operator
// *************************************************************

template <class T>
Histogram<T> &	Histogram<T>::operator= (const Histogram<T> &other)
{
	if (this==&other)
		return *this;

	// Sizes do not match -> realloc
	if (_size != other._size)
	{
		delete bins;
		_size = other._size;
		bins = new T [_size];
	}

	for (int i=0; i<_size; ++i)
		bins[i] = other.bins[i];

	return *this;
}

// *********************************************************************
// Clear the Histogram
// *********************************************************************

template <class T>
inline void Histogram<T>::clear ()
{
	for	(int i=0; i<_size; ++i)
		bins[i]=0;
}

// *********************************************************************
// Get the histogram bin index from the value
// *********************************************************************

template <class T>
int Histogram<T>::value2index (T val)
{
	int index = (int) ((double) (val - min) / step) / 1;
	if (index>=_size) index = _size-1;
	if (index<0) index = 0;
	return index;
}

// *********************************************************************
// Add a value
// *********************************************************************

template <class T>
inline void Histogram<T>::add (T val)
{
	++bins[value2index(val)];
}

// *********************************************************************
// Add a value - by weight
// *********************************************************************

template <class T>
inline void Histogram<T>::add (T val, T weight)
{
	bins[value2index(val)] += weight;
}


// *********************************************************************
// Remove a value
// *********************************************************************

template <class T>
void Histogram<T>::remove (T val)
{
	--(bins[value2index(val)]);
}

// ******************************************************************
// Gets	the	sum	of the histogram
// ******************************************************************

template <class T>
inline T Histogram<T>::sum () const
{
	T rv=0;
	for	(int i=0; i<_size; ++i)
		rv += bins[i];
	return rv;
}

// ******************************************************************
// Gets	the	mean
// ******************************************************************

template <class T>
double Histogram<T>::getMean()
{
	double total=0;
	double run_val=min;
	for (int i=0; i<_size; ++i)
	{
		total += bins[i]*run_val;
		run_val += step;
	}
	return sum() == 0 ? 0 : total / (double) sum();
}

// ******************************************************************
// Gets	the variance
// ******************************************************************

template <class T>
double Histogram<T>::getVariance()
{
	double v, total=0, run_val=min, m=getMean();
	for (int i=0; i<_size; ++i)
	{
		v=(run_val-m);
		total += bins[i]*v*v;
		run_val += step;
	}
	return sum() == 0 ? 0 : total / (double) sum();
}

// *********************************************************************
// Normalize the histogram
// *********************************************************************

template <class T>
inline void Histogram<T>::normalize ()
{
	T masse = sum();
	for	(int i=0; i<_size; ++i)
	{
		bins[i]	/= masse;
	}
}

// *********************************************************************
// - operator: L2 norm
// *********************************************************************

template <class T>
inline double Histogram<T>::operator - (const Histogram<T> & other)
{
 	double rv=0;
	for (int i=0; i<_size; ++i)
		rv += ((bins[i]-other.bins[i]) * (bins[i]-other.bins[i]));
	return rv;
}

// *********************************************************************
// The difference histogram
// *********************************************************************

template <class T>
Histogram<T> * Histogram<T>::diffHist (const Histogram<T> &other)
{
	Histogram<T> *rv = new Histogram<T> (*this);
	for (int i=0; i<_size; ++i)
		rv->bins[i] -= other.bins[i];
	return rv;
}
// *********************************************************************
// Get the index of the maximum
// *********************************************************************

template <class T>
int Histogram<T>::argmax() {
	bool first=true;
	int argmax=0;
	T max=0;
	for (int i=0; i<_size; ++i) {
		if (first || bins[i] > max) {
			max = bins[i];
			argmax = i;
			first = false;
		}
	}
	return argmax;
}

// *********************************************************************
// Get the value of the bin center
// *********************************************************************

template <class T>
inline T Histogram<T>::getBinCenter(int bin) {
	return (T) min + step * ((double) bin + 0.5);
}

// *********************************************************************
// Get borders
// *********************************************************************

template <class T>
inline int Histogram<T>::firstNonZero() {
   for (int i=0; i<_size; ++i)
   	if (bins[i]!=0) return i;
   return _size-1;
}

template <class T>
inline int Histogram<T>::lastNonZero() {
   for (int i=_size-1; i>0; --i)
   	if (bins[i]!=0) return i;
   return 0;
}


// *********************************************************************
// Get the index of percentage of pixels left
// *********************************************************************

template <class T>
inline int Histogram<T>::getIndexPercLeft (double perc) {

	int i;
	T M, si, border_mass;

	M = sum();

	if (M==0)
		return 0;

	border_mass = (int) ((perc * M) / 100);

	if (border_mass == 0)
		border_mass = 1;

	// Search left border
	si = 0;
	for	(i=0; i<_size; ++i) {
		if (si >= border_mass)
			break;
		si += bins[i];
	}

	return i;
}

// *********************************************************************
// Get the index of percentage of pixels right
// *********************************************************************


template <class T>
inline int Histogram<T>::getIndexPercRight (double perc) {
	int i;
	T si;
	T M, border_mass;

	M = sum();

	if (M==0)
		return _size-1;

	border_mass = (T)((perc * M) / 100.0);

	if (border_mass == 0)
		border_mass = 1;

	// Search right border
	si = 0;
	for	(i=_size-1; i>=0; --i) {
		if (si >= border_mass)
			break;
		si += bins[i];
	}

	return i;
}

// *********************************************************************
// Cut an interval left and right

// *********************************************************************

template <class T>
inline void Histogram<T>::cutOutliers (double leftperc, double rightperc) {
	int left_border, right_border;

	left_border = getIndexPercLeft (leftperc);
	right_border = getIndexPercRight (rightperc);

	// Truncate the histogram
	for (int i=0; i<= left_border && i<_size; ++i) bins[i] = 0;
	for (int i=_size-1; i>=right_border && i>0; --i) bins[i] = 0;
}

// *********************************************************************
// A special feature for histograms which store grayvalues:
// Get the local contrast, i.e. the span of the gray values
// Get the minimum and maximum of the histogram,
// But remove the first and the last 12% of the
// values to make the algo more robust to outliers
// *********************************************************************

template <class T>
inline unsigned char Histogram<T>::getContrast () {
	int left_border, right_border, contrast;
	int si, border_mass;
	int M, i;

	M = sum();
	border_mass = (15 * M) / 100;

	if (border_mass == 0)
		border_mass = 1;

	// Search left border
	si = 0;
	for	(i=0; i<_size; ++i) {
		if (si >= border_mass)
			break;
		si += bins[i];
	}
	left_border = i;

	// Search right border
	si = 0;
	for	(i=_size-1; i>=0; --i) {
		if (si >= border_mass)
			break;
		si += bins[i];
	}
	right_border = i;

	contrast = right_border - left_border;
	if (contrast<0)
		contrast = 0;

	return contrast;
}

// *************************************************************
// Compute the optimal threshold value with the fisher method.
// The Implementation gathered from Jean-Michel Jolion
// The last two arguments are filled with mode peaks, if they
// are not null (the arguments).
// *************************************************************

template <class T>
int Histogram<T>::getThresholdValueFisher (int kmin, int kmax,
	float *out_peakleft, float *out_peakright) {

	int	k, kopt ;
	float W, m0, m1, M, S, Smax ;
	float m1opt, m0opt;
	T *H = bins;

	if (kmin<0) kmin = 0;
	if (kmax>255) kmax = 255;

	// Search the borders
	kopt = 0 ;
	m1 = W = 0.	;
	if (kmin ==	0)
		while (H[kmin] == 0)
			kmin++ ;
	if (kmax >=	_size-1)
		while (H[kmax]	== 0)
			kmax-- ;
	else {
	 	if	(kmax == -1) {
	 		kmax = kmin ;
	 		Smax =	H[kmin]	;

	   		for (k = kmin+1 ; k < _size-1 ; k++) {

				if (H[k] > Smax) {
					Smax	= H[k] ;
					kmax =	k ;
				}
			}
		}
	}
	while (H[kmin] == 0)
		kmin++	;
	while (H[kmax] == 0)
		kmax--	;

	// Masse -> M
	M =	0 ;
	for	(k = kmin ;	k <= kmax ;	k++)
		M += H[k] ;

	m0 = m0opt = kmin ;
	for	(k = kmin+1	; k	<= kmax	; k++) {
		m1 +=	k*H[k] ;
	  	W	+= H[k]	;
	}

	m1 = m1opt = m1/W ;
	kopt = kmin	;
	Smax = (W/M)*(1.-W/M)*(m1 -	m0)*(m1	- m0) ;

	for	(k = kmin+1	; k	< kmax ; k++) {
		m0 = m0*(M - W) +	k*H[k] ;
	  	m1 = m1 *	W -	k*H[k] ;
	  	W	-= H[k]	;
	  	m1 /=	W ;
	  	m0 /=	(M-W) ;
	  	S	= (W/M)*(1.-W/M)*(m1 - m0)*(m1 - m0) ;

	  	if (S	> Smax)	{
	  		Smax = S ;
	  		kopt =	k ;
	  		m0opt = m0;
	  		m1opt = m1;
	  	}
	}

#ifdef DEBUG
		/*
		cerr << "threshold by Fisher: " << kopt << endl;
		cerr << "m1 = " << m1opt << ", m0 = " << m0opt << endl;
		*/
#endif

	if (out_peakleft != NULL) *out_peakleft = m0;
	if (out_peakright != NULL) *out_peakright = m1;
	return kopt;
}

// *************************************************************
// Compute the optimal threshold value with the fisher method.
// The TWO THRESHOLD VERSION!!!!
// *************************************************************

#define FISH_CRIT	((W0/M)*(W1/M)*(m1-m0)*(m1-m0)+(W1/M)*(W2/M)*(m2-m1)*(m2-m1))
#define FISH_OUT	; // {	if (1) cerr << k1 << " " << k2 << " " << S << endl; }

template <class T>
void Histogram<T>::getThresholdValueFisher2th (int kmin, int kmax, int &out_k1, int &out_k2)
{
	int	k1, k2, k1opt, k2opt;
	float W0, W1, W2, m0, m1, m2, M, S, Smax ;
	float m0opt, m1opt, m2opt;
	T *H = bins;

	// Search the borders. Modes with zero contents produce divisions by zero.
	if (kmin ==	0)
		while (H[kmin] == 0)
			kmin++ ;
	if (kmax >=	_size-1)
	{
		kmax = _size-1;
		while (H[kmax]	== 0)
			kmax-- ;
	}

	while (H[kmin] == 0)
		kmin++	;
	while (H[kmax] == 0)
		kmax--	;

	// Masse -> M
	M =	0 ;
	for	(int k = kmin ;	k <= kmax ;	k++)
		M += H[k] ;

	// -------------------------------------------------------------
	// Initialization:
	// -------------------------------------------------------------

	// Mode 0:
	m0 = m0opt = kmin ;
	W0 = H[kmin];

	// Mode 1:
	m1 = m1opt = kmin+1;
	W1 = H[kmin+1];

	// Mode 2:
	for	(k2 = kmin+2; k2 <= kmax; k2++)
	{
		m2 += k2*H[k2] ;
	  	W2 += H[k2]	;
	}
	m2 = m2opt = m2/W2 ;

	// Thresholds and criterion
	k1opt = k1 = kmin;
	k2opt = k2 = kmin+1;
	S = Smax = FISH_CRIT;

	FISH_OUT;

	// -------------------------------------------------------------
	// Travers the search space
	// Move the first threshold
	// -------------------------------------------------------------
	k1 = kmin+1;
	while (k1<kmax-1) {

		// Recalculate the statistics of mode 0 from the last round
		m0 = m0*W0 + k1*H[k1] ;
		W0 = W0 + H[k1];

		// W0 can never be zero, because at least the first bin is
		// non-zero
		m0 /= W0;

		// Search the first non-zero bin, a zero mode is useless and leads to
		// division-by-zeros. There must be a non-zero bin, since we set the
		// right border at the beginning of the algorithm.
		k2 = k1+1;
		while (k2<kmax && H[k2]==0)
			++k2;
		if (k2==kmax)
		{
			++k1;
			continue;
		}

		// Initialize modes nr. 1 and 2 for the inner loop (second threshold)
		m1=k2;
		W1=H[k2];
		m2=W2=0;
		for	(int i = k2+1; i<= kmax; ++i)
		{

			m2 += i*H[i] ;
	  		W2 += H[i]	;
		}
		m2 /= W2;

		S = FISH_CRIT;
		if (S>Smax)
		{
    	  		Smax = S ;
    	  		k1opt =	k1;
    	  		k2opt = k2;
    	  		m0opt = m0;
    	  		m1opt = m1;
    	  		m2opt = m2;
    	}

    	FISH_OUT;

    	// -------------------------------------------------------------
    	// Move the second threshold
    	// -------------------------------------------------------------

		++k2;
		while (k2<kmax)
		{
    		m1 = m1*W1 + k2*H[k2] ;
    	  	m2 = m2*W2 - k2*H[k2] ;
    	  	W1	+= H[k2]	;
    	  	W2	-= H[k2]	;
    	  	m1 /= W1 ;
    	  	m2 /= W2;

    	  	S = FISH_CRIT;
    		if (S>Smax)
    		{
        	  		Smax = S ;
        	  		k1opt =	k1;
        	  		k2opt = k2;
        	  		m0opt = m0;
        	  		m1opt = m1;
        	  		m2opt = m2;
        	}

        	FISH_OUT;

    	  	++k2;
		}
		++k1;
	}

	cerr << "thresholds by 2-th Fisher: " << k1opt << "," << k2opt << endl;
	cerr << "m0 = " << m0opt << ", m1 = " << m0opt
		 << ", m2 = " << m2opt << "\n" << endl;

	out_k1 = k1opt;
	out_k2 = k2opt;
}


// *********************************************************************
// Get the intraclass variance for all pixels of both classes
// separated by threshold t
// *********************************************************************

template <class T>
double Histogram<T>::getIntraclassVariance (int t) {
	double mass0=0,
		   mass1=0;
	double mean0=0,
		   mean1=0;
	double x, var=0;

	// Calculate the mean of the two classes
	for (int i=0; i<t; ++i) {

			mass0 += bins[i];

			mean0 += i*bins[i];
	}
	for (int i=t; i<_size; ++i) {
			mass1 += bins[i];
			mean1 += i*bins[i];
	}
	mean0 /= mass0;
	mean1 /= mass1;

	for (int i=0; i<t; ++i) {
		x = mean0 - bins[i];
		var += x*x;
	}
	for (int i=t; i<_size; ++i) {
		x = mean1 - bins[i];
		var += x*x;
	}
	var /= (mass0+mass1);

	return var;
}

// *********************************************************************
// Perform a histogram stretch
// *********************************************************************

template <class T>
inline void Histogram<T>::stretch() {

	int omin, omax;
	int nmin, nmax;

	omin = firstNonZero();
	omax = lastNonZero();
	nmin = 0;
	nmax = _size-1;

	/*

	// Travers all bins
	for (int i=min; i<=max; ++i) {

    	val -= omin;
    	val = (val*(nmax-nmin))/(omax-omin);
    	val += nmin;

    	if (val>nmax)
    		val = nmax;
    	else
    		if (val < nmin)
    			val = nmin;
	}

	*/
}

// *********************************************************************
// Remove false accumulations in bin 0 and _size-1
// *********************************************************************

template <class T>
void Histogram<T>::removeFalsePeaks() {

	bins[0] = bins[1];
	bins[_size-1] = bins[_size-2];

	/*
	double avg;

	avg=0;
	for (i=0; i<5; ++i)
		avg += bins[i];
	avg /= 5.0;

	if (((double) bins[0])/avg > 2.0)
		bins[0] = bins[1];

	avg=0;
	for (i=0; i<5; ++i)
		avg += bins[_size-1-i];
	avg /= 5.0;

	if (((double) bins[_size-1])/avg > 2.0)
		bins[_size-1] = bins[_size-2];
	*/
}

// *********************************************************************
// Transform the histogram into a cumulative histogram
// *********************************************************************

template<class T>
void Histogram<T>::toCumulative()
{
	T sum=0;
	for (int i=0; i<_size; ++i)
	{
		sum+=bins[i];
		bins[i] = sum;
	}
}

template<class T>
void Histogram<T>::toCumulativeFromRight()
{
	T sum=0;
	for (int i=_size-1; i>-1; --i)
	{
		sum+=bins[i];
		bins[i] = sum;
	}
}


template<class T>
void Histogram<T>::toLogScale()
{
	for (int i=0; i<_size; ++i)
	{
	    if(bins[i]>0)
            bins[i] = (T)log(double(bins[i]));
	}
}

template<class T>
void Histogram<T>::truncFromRightStop(double percentage)
{
    assert(percentage <=1.0);

    normalize(); // to make sure the sum=1

	double sum=0;
	for (int i=_size-1; i>-1; --i)
	{
		if(sum>percentage)
            bins[i] = (T)0;

        sum+=bins[i];
	}
}
// *********************************************************************
// Perform the Kolmogorov-Smirnoff statistical test.
// ATTENTION:
// The two histograms MUST BE CUMULATIVE HISTOGRAMS!!!!!!!!!
// *********************************************************************

template<class T>
float Histogram<T>::kolmogorovSmirnovTest(Histogram<T> &other) const
{
	T D,Dmax=0;

	for (int i=0; i<_size; ++i)
	{
    	D=fabs(bins[i] - other.bins[i]) ;
    	if (D > Dmax)
        	Dmax = D;
	}

	return Dmax;
}

// *********************************************************************
// Create a uniform statistical distribution
// *********************************************************************

template<class T>
void Histogram<T>::createUniformDistributionFromValue (T val)
{
	for (int i=0; i<_size; ++i)
		bins[i] = val;
}

template<class T>
void Histogram<T>::createUniformDistributionFromMass (T mass)
{
	sum=0;
	for (int i=0; i<_size; ++i)
		bins[i] += (mass*(i+1))/_size-sum;
}

// *************************************************************
// Filter the data with a given filter mask.
// The mask must of course be 1D
// *************************************************************
/*
template<class T>
void Histogram<T>::filter (DataSerie &fm)
{
	T *FIL, *todelete;
	int	ix, bx, fx;
	double newVal;
	int fmx = fm.size();

	// Make a copy of the data
	FIL = new T[_size];

	// Treat the borders
	for (int i=0; i<fmx/2; ++i)
	{
		FIL[i] = bins[i];
		FIL[_size-1-i] = bins[_size-1-i];
	}

	// Travers the data
	for	(int x=0+fmx/2; x<(_size-fmx/2); ++x )
	{
		// the left border
		bx = x-(fmx/2);

    	// convolve
    	newVal = 0;
        for	(ix=bx,	fx=0; fx<fmx; ++ix, ++fx)
			newVal += (((double)bins[ix])*fm[fx]);

		FIL[x]=	(int) rint (newVal);
	}

	// Copy result and clean up
	todelete = bins;
	bins = FIL;
	delete todelete;
}
*/
// *********************************************************************
// OUTPUT - normal
// *********************************************************************

template <class T>
inline ostream & operator << (ostream & os, Histogram<T> h)
{
	double class_med;

	class_med = h.getMin();
	if (!h.discret)
		class_med += 0.5*h.getStepSize();
	for (int i=0; i<h._size; ++i) {
		os << class_med << " " << h[i] << endl;
		class_med += h.getStepSize();
	}
	return os;
}

// A special case for double: the compiler can't seem to find it :(
/*
inline ostream & operator << (ostream & os, Histogram<double> h)
{
	double class_med;

	class_med = h.min;
	if (!h.discret)
		class_med += 0.5*h.step;
	for (int i=0; i<h._size; ++i) {
		os << class_med << " " << h[i] << endl;
		class_med += h.step;
	}
	return os;
}

// A special case for byte: needs type casting

template <byte>
ostream & operator << (ostream & os, Histogram<byte> h) {
	for (int i=h.firstNonZero(); i<=h.lastNonZero(); ++i)
		os << i << " " << (int) h[i] << endl;
	return os;
}
*/

// A special case for int: needs type casting
inline ostream & operator << (ostream & os, Histogram<int> &h)
{
	for (int i=h.firstNonZero(); i<=h.lastNonZero(); ++i)
		os << i << " " << (int) h[i] << endl;
	return os;
}

// *********************************************************************
// OUTPUT - for matlab
// *********************************************************************

template <class T>
inline void Histogram<T>::matlabOut (ostream & os, const char *s)
{
	os << s << " = [";
	for (int i=0; i<_size; ++i)
		os << " " << (T) bins[i];
	os << "];\n";
}

// *********************************************************************
// OUTPUT - for gnuplot
// *********************************************************************

template <class T>
inline void Histogram<T>::gnuplotOut (  const char *title,
                                        const char *Xlabel,
                                        const char *Ylabel,
                                        bool log_scale)
{
    FILE * pFile = fopen("./res_gnuplot/histogram.txt","w"); // empty results file
    if (!pFile)
    {
        cout << "Cannot open the file: ./res_gnuplot/histogram.txt!!!" << endl;
        return;
    }

    normalize ();

    double xrangemin=1e16, xrangemax=-1e16;
    char sentence [256];
	int i=0;
	for (; i < _size; ++i)
	{
	    if(getBinCenter(i)<xrangemin) xrangemin = getBinCenter(i);

	    if(getBinCenter(i)>xrangemax) xrangemax = getBinCenter(i);

	    sprintf (sentence, "%f %f\n", getBinCenter(i), (T) bins[i]);
        fputs (sentence,pFile);
	}

    //sprintf (sentence, "%d %f\n", getBinCenter(i), (T) 0);
    //fputs (sentence,pFile);

    //cout << "xrangemin = " << xrangemin << " xrangemax = " << xrangemax << endl;

    fclose (pFile);

    gnuplot_interface Histo("./res_gnuplot/histogram.txt", title, Xlabel, Ylabel,"lines");
    Histo.gnuplot_histo_png((char*)"./res_gnuplot/images/resHisto.png", xrangemin, xrangemax, log_scale);
}
#endif

