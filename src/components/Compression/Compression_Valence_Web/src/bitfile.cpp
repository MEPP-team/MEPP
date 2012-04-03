/***************************************************************************
*                       Bit Stream File Class Implementation
*
*   File    : bitfile.cpp
*   Purpose : This file implements a simple class of I/O methods for files
*             that contain data in sizes that aren't integral bytes.  An
*             attempt was made to make the methods in this class analogous
*             to the methods provided to manipulate file streams.  The
*             methods contained in this class were created with compression
*             algorithms in mind, but may be suited to other applications.

****************************************************************************

/***************************************************************************
*                             INCLUDED FILES
***************************************************************************/
#include "bitfile.h"

using namespace std;

/***************************************************************************
*                            TYPE DEFINITIONS
***************************************************************************/

/* union used to test for endianess */
typedef union
{
    unsigned long word;
    unsigned char bytes[sizeof(unsigned long)];
} endian_test_t;

/***************************************************************************
*                                 METHODS
***************************************************************************/

/***************************************************************************
*   Method     : bit_file_c - default constructor
*   Description: This is the default bit_file_c constructor.  It
*                initializes stream pointers to NULL and clears the bit
*                buffer.
*   Parameters : None
*   Effects    : Initializes private members.
*   Returned   : None
***************************************************************************/
bit_file_c::bit_file_c(void)
{
    m_InStream = NULL;
    m_OutStream = NULL;
    m_BitBuffer = 0;
	m_UnsignedBuffer=0;
    m_BitCount = 0;
    m_Mode = BF_NO_MODE;
	m_Current_Decoded_Bits = 0;

    /* test for endianess */
    endian_test_t endianTest;

    endianTest.word = 1;

    if (endianTest.bytes[0] == 1)
    {
        /* LSB is 1st byte (little endian)*/
        m_endian = BF_LITTLE_ENDIAN;
    }
    else if (endianTest.bytes[sizeof(unsigned long) - 1] == 1)
    {
        /* LSB is last byte (big endian)*/
        m_endian = BF_BIG_ENDIAN;
    }
    else
    {
        m_endian = BF_UNKNOWN_ENDIAN;
    }
}

/***************************************************************************
*   Method     : bit_file_c - constructor
*   Description: This is a bit_file_c constructor.  It opens an input or
*                output stream and clears the bit buffer.  An exception
*                will be thrown on error.
*   Parameters : fileName - NULL terminated string containing the name of
*                           the file to be opened.
*                mode - The mode of the file to be opened
*   Effects    : Initializes private members.  Creates and opens an input
*                or output stream.
*   Returned   : None
*   Exception  : "Error: Invalid File Type" - for unknown mode
*                "Error: Unable To Open File" - if stream cannot be opened
***************************************************************************/
bit_file_c::bit_file_c(const char *fileName, const BF_MODES mode)
{
    m_InStream = NULL;
    m_OutStream = NULL;
    m_BitBuffer = 0;
    m_BitCount = 0;

    switch (mode)
    {
        case BF_READ:
            m_InStream = new ifstream(fileName, ios::in | ios::binary);

            if (!m_InStream->good())
            {
                delete m_InStream;
                m_InStream = NULL;
            }
            else
            {
                m_Mode = mode;
            }
            break;

        case BF_WRITE:
            m_OutStream = new ofstream(fileName, ios::out | ios::binary);

            if (!m_OutStream->good())
            {
                delete m_OutStream;
                m_OutStream = NULL;
            }
            else
            {
                m_Mode = mode;
            }
            break;

        case BF_APPEND:
            m_OutStream =
                new ofstream(fileName, ios::out | ios::binary | ios::app);

            if (!m_OutStream->good())
            {
                delete m_OutStream;
                m_OutStream = NULL;
            }
            else
            {
                m_Mode = mode;
            }
            break;

        default:
            throw("Error: Invalid File Type");
            break;
    }

    /* make sure we opened a file */
    if ((m_InStream == NULL) && (m_OutStream == NULL))
    {
        throw("Error: Unable To Open File");
    }

    /* test for endianess */
    endian_test_t endianTest;

    endianTest.word = 1;

    if (endianTest.bytes[0] == 1)
    {
        /* LSB is 1st byte (little endian)*/
        m_endian = BF_LITTLE_ENDIAN;
    }
    else if (endianTest.bytes[sizeof(unsigned long) - 1] == 1)
    {
        /* LSB is last byte (big endian)*/
        m_endian = BF_BIG_ENDIAN;
    }
    else
    {
        m_endian = BF_UNKNOWN_ENDIAN;
    }
}

/***************************************************************************
*   Method     : ~bit_file_c - destructor
*   Description: This is the bit_file_c destructor.  It closes and frees
*                any open file streams.  The bit buffer will be flushed
*                prior to closing an output stream.
*   Parameters : None
*   Effects    : Closes and frees open file streams.
*   Returned   : None
***************************************************************************/
bit_file_c::~bit_file_c(void)
{
    if (m_InStream != NULL)
    {
        m_InStream->close();
        delete m_InStream;
    }

    if (m_OutStream != NULL)
    {
        /* write out any unwritten bits */
        if (m_BitCount != 0)
        {
            m_BitBuffer <<= (8 - m_BitCount);
            m_OutStream->put(m_BitBuffer);
        }

        m_OutStream->close();
        delete m_OutStream;
    }
}

/***************************************************************************
*   Method     : Open
*   Description: This method opens an input or output stream and initializes
*                the bit buffer.  An exception will be thrown on error.
*   Parameters : fileName - NULL terminated string containing the name of
*                           the file to be opened.
*                mode - The mode of the file to be opened
*   Effects    : Creates and opens an input or output stream.
*   Returned   : None
*   Exception  : "Error: File Already Open" - if object has an open file
*                "Error: Invalid File Type" - for unknown mode
*                "Error: Unable To Open File" - if stream cannot be opened
***************************************************************************/
void bit_file_c::Open(const char *fileName, const BF_MODES mode)
{
    /* make sure file isn't already open */
    if ((m_InStream != NULL) || (m_OutStream != NULL))
    {
        throw("Error: File Already Open");
    }

    switch (mode)
    {
        case BF_READ:
            m_InStream = new ifstream(fileName, ios::in | ios::binary);

            if (!m_InStream->good())
            {
                delete m_InStream;
                m_InStream = NULL;
            }
            else
            {
                m_Mode = mode;
            }

            m_BitBuffer = 0;
            m_BitCount = 0;
            break;

        case BF_WRITE:
            m_OutStream = new ofstream(fileName, ios::out | ios::binary);

            if (!m_OutStream->good())
            {
                delete m_OutStream;
                m_OutStream = NULL;
            }
            else
            {
                m_Mode = mode;
            }

            m_BitBuffer = 0;
            m_BitCount = 0;
            break;

        case BF_APPEND:
            m_OutStream =
                new ofstream(fileName, ios::out | ios::binary | ios::app);

            if (!m_OutStream->good())
            {
                delete m_OutStream;
                m_OutStream = NULL;
            }
            else
            {
                m_Mode = mode;
            }

            m_BitBuffer = 0;
            m_BitCount = 0;
            break;

        default:
            throw("Error: Invalid File Type");
            break;
    }

    /* make sure we opened a file */
    if ((m_InStream == NULL) && (m_OutStream == NULL))
    {
        throw("Error: Unable To Open File");
    }
}

/***************************************************************************
*   Method     : Close
*   Description: This method closes and frees any open file streams.  The
*                bit buffer will be flushed prior to closing an output
*                stream.  All member variables are re-initialized.
*   Parameters : None
*   Effects    : Closes and frees open file streams.  Resets member
*                variables.
*   Returned   : None
***************************************************************************/
void bit_file_c::Close(void)
{
    if (m_InStream != NULL)
    {
        m_InStream->close();
        delete m_InStream;

        m_InStream = NULL;
        m_BitBuffer = 0;
		m_UnsignedBuffer = 0;
        m_BitCount = 0;
        m_Mode = BF_NO_MODE;
    }

    if (m_OutStream != NULL)
    {
        /* write out any unwritten bits */
        if (m_BitCount != 0)
        {
            m_UnsignedBuffer <<= (8 - m_BitCount);
            m_OutStream->put(m_UnsignedBuffer);
        }

        m_OutStream->close();
        delete m_OutStream;

        m_OutStream = NULL;
        m_BitBuffer = 0;
        m_BitCount = 0;
        m_Mode = BF_NO_MODE;
    }
}





///***************************************************************************
//*   Method     : decode_bits
//*   Description: This method reads the specified number of bits from the
//*                input stream and writes them to the requested memory
//*                location (msb to lsb).
//*   Parameters : 
//*                count - number of bits to read
//*   Effects    : Reads bits from the bit buffer and file stream.  The bit
//*                buffer will be modified as necessary.
//*   Returned   : EOF for failure, otherwise the number of bits read.  If
//*                an EOF is reached before all the bits are read, bits
//*                will contain every bit through the last complete byte.
//***************************************************************************/

int bit_file_c::decode_bits(const unsigned int count)
{
	//update the total number of bits decoded
	m_Current_Decoded_Bits += count;
	
	//decode the value
    unsigned long returnValue=0;
	while(count>m_BitCount)
	{m_UnsignedBuffer<<= 8 ;
	m_UnsignedBuffer|=m_InStream->get();
	m_BitCount +=8;
	}
	returnValue = m_UnsignedBuffer >> (m_BitCount-count);

	m_UnsignedBuffer -= (returnValue << (m_BitCount-count));
	m_BitCount-=count;
	return returnValue;

}
///***************************************************************************
//*   Method     : encode_bits
//*   Description: This method writes the bit passed as a parameter to the
//*                output stream.
//*   Parameters : bits - the bit value to be written
//*				   count - number of bits to be written 
//*   Effects    : Writes a bit to the bit buffer.  If the buffer has a byte,
//*                the buffer is written to the output stream and cleared.
//*   Returned   : On success, the bit value written, otherwise EOF.
//***************************************************************************/

int bit_file_c::encode_bits(int bits, const unsigned int count)
{
	//move the buffer count times
	m_UnsignedBuffer <<= (count);
	//xor the new arriving bits to the buffer
	unsigned mask = 0xffffffff;
	m_UnsignedBuffer|= (bits&(mask>>(32-count)));
	//update the total number of bits to be encoded
	m_BitCount+=count;

	int value_to_write;
	//write the bytes in buffer 
	while(m_BitCount>=8)
	{
		value_to_write = m_UnsignedBuffer>>(m_BitCount-8);
		m_OutStream->put((char)value_to_write);

		//update the buffer
		m_UnsignedBuffer -= (value_to_write<<(m_BitCount-8));
		//update the number of remaining bits
		m_BitCount-=8;
       
	}

	return 1;


}

unsigned int bit_file_c::calculate_current_decoded_size()
{
	return m_Current_Decoded_Bits;
}

bool bit_file_c::eof(void)
{
    if (m_InStream != NULL)
    {
        return (m_InStream->eof());
    }

    if (m_OutStream != NULL)
    {
        return (m_OutStream->eof());
    }

    /* return false for no file */
    return false;
}

/***************************************************************************
*   Method     : good
*   Description: This method is analogous to good for file streams.
*   Parameters : None
*   Effects    : None
*   Returned   : Returns good for the opened file stream.  False is
*                returned if there is no open file stream.
***************************************************************************/
bool bit_file_c::good(void)
{
    if (m_InStream != NULL)
    {
        return (m_InStream->good());
    }

    if (m_OutStream != NULL)
    {
        return (m_OutStream->good());
    }

    /* return false for no file */
    return false;
}

/***************************************************************************
*   Method     : bad
*   Description: This method is analogous to bad for file streams.
*   Parameters : None
*   Effects    : None
*   Returned   : Returns bad for the opened file stream.  False is
*                returned if there is no open file stream.
***************************************************************************/
bool bit_file_c::bad(void)
{
    if (m_InStream != NULL)
    {
        return (m_InStream->bad());
    }

    if (m_OutStream != NULL)
    {
        return (m_OutStream->bad());
    }

    /* return false for no file */
    return false;
}
