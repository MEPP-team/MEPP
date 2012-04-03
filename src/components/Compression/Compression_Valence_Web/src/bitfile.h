/***************************************************************************
*                       Bit Stream File Class Header
*
*   File    : bitfile.h
*   Purpose : Provides definitions and prototypes for a simple class of I/O
*             methods for files that contain data in sizes that aren't
*             integral bytes.  An attempt was made to make the methods in
*             this class analogous to the methods provided to manipulate
*             file streams.  The methods contained in this class were
*             created with compression algorithms in mind, but may be
*             suited to other applications.
*   Author  : Michael Dipperstein
*   Date    : July 20, 2004
*
****************************************************************************
*   UPDATES
*
*   $Id: bitfile.h,v 1.6 2008/01/27 06:04:54 michael Exp $
*   $Log: bitfile.h,v $
*   Revision 1.6  2008/01/27 06:04:54  michael
*   Added  ByteAlign() and FlushOutput() methods.
*
*   Revision 1.5  2007/08/26 21:41:36  michael
*   All methods that don't modify the calling object have been made const
*   to increase functionality of const bit_array_c.
*
*   Changes required for LGPL v3.
*
*   Revision 1.4  2007/07/16 02:07:16  michael
*   Use -pedantic option when compiling.
*
*   Revision 1.3  2005/12/10 05:20:01  michael
*   Added methods to get/put bits from/to integer types.
*
*   Revision 1.2  2005/06/23 04:39:06  michael
*   Convert from DOS end of line to Unix end of line
*
*   Revision 1.1.1.1  2004/08/04 13:45:38  michael
*   bitfile class
*
*
****************************************************************************
*
* Bitfile: Bit Stream File I/O Class

***************************************************************************/

#ifndef __BITFILE_H
#define __BITFILE_H

#include <iostream>
#include <fstream>
#include <math.h>

/***************************************************************************
*                            TYPE DEFINITIONS
***************************************************************************/
typedef enum
{
    BF_READ = 0,
    BF_WRITE = 1,
    BF_APPEND= 2,
    BF_NO_MODE
} BF_MODES;

typedef enum
{
    BF_UNKNOWN_ENDIAN,
    BF_LITTLE_ENDIAN,
    BF_BIG_ENDIAN
} endian_t;

class bit_file_c
{
    public:
        bit_file_c(void);
        bit_file_c(const char *fileName, const BF_MODES mode);
        virtual ~bit_file_c(void);

        /* open/close bit file */
        void Open(const char *fileName, const BF_MODES mode);
        void Close(void);

 
        /* fill byte with ones or zeros and write out results */
 


 		int decode_bits( const unsigned int count);
        int encode_bits(int bits, const unsigned int count);
		//int encode_bits_2(int bits, const unsigned int count);
		unsigned int calculate_current_decoded_size();
        /* status */
        bool eof(void);
        bool good(void);
        bool bad(void);

    private:
        std::ifstream *m_InStream;      /* input file stream pointer */
        std::ofstream *m_OutStream;     /* output file stream pointer */
        endian_t m_endian;              /* endianess of architecture */
        char m_BitBuffer;               /* bits waiting to be read/written */
        unsigned int m_UnsignedBuffer;               /* bits waiting to be read/written */
        unsigned char m_BitCount;       /* number of bits in bitBuffer */
        BF_MODES m_Mode;                /* open for read, write, or append */
		unsigned int m_Current_Decoded_Bits;       /* number of decoded bits in total */
        

       
};

#endif  /* ndef __BITFILE_H */
