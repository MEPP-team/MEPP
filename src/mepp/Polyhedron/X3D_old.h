/*
 * X3DMeshExtractor extract IndexedFaceSet from X3D document
 * and store it into IFSData object in memory and/or hard drive.
 * IFSData object includes faces, vertices, colors and infos. 
 */

#ifndef X3DMESHEXTRACTOR_H_
#define X3DMESHEXTRACTOR_H_

// INCLUDES //

// STL includes
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

// libc includes
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <fcntl.h>

// Main includes
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/TransService.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>

// FilterHandlers includes
#include <xercesc/parsers/SAX2XMLFilterImpl.hpp>
#include <xercesc/sax2/Attributes.hpp>

// Handlers includes
#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/framework/XMLFormatter.hpp>

// Personnal includes
#include "X3D_IFS_old.h"

// MT
#ifdef _MSC_VER
#pragma warning(disable : 4355)
#endif

XERCES_CPP_NAMESPACE_USE

// DEFINITIONS //

// Max size of a line read from a file
#define BUFFER_SIZE 256

// State machine values
#define NONE           0
#define SHAPE          1
#define INDEXEDFACESET 2
#define COORDINATE     3
#define COLOR          4

// STATIC CONST DATA //

static const XMLCh gEndElement[] =
{	chOpenAngle, chForwardSlash, chNull};

static const XMLCh gEndPI[] = { chQuestion, chCloseAngle, chNull };

static const XMLCh gStartPI[] = { chOpenAngle, chQuestion, chNull };

static const XMLCh gXMLDecl1[] = { chOpenAngle, chQuestion, chLatin_x,
		chLatin_m, chLatin_l, chSpace, chLatin_v, chLatin_e, chLatin_r,
		chLatin_s, chLatin_i, chLatin_o, chLatin_n, chEqual, chDoubleQuote,
		chDigit_1, chPeriod, chDigit_0, chDoubleQuote, chSpace, chLatin_e,
		chLatin_n, chLatin_c, chLatin_o, chLatin_d, chLatin_i, chLatin_n,
		chLatin_g, chEqual, chDoubleQuote, chNull };

static const XMLCh gXMLDecl2[] = { chDoubleQuote, chQuestion, chCloseAngle,
		chLF, chNull };

// CLASSES //

struct Attr {
	const XMLCh* qName;
	const XMLCh* uri;
	const XMLCh* localPart;
	const XMLCh* value;
	const XMLCh* attrType;
};

class AttrList : public Attributes, public RefVectorOf<Attr> {

public:

	AttrList(unsigned count) :
		RefVectorOf<Attr>(count) {
	}

	virtual /*unsigned int*/XMLSize_t getLength() const {
		return size();
	}

	virtual const XMLCh* getURI(const /*unsigned int*/XMLSize_t index) const {
		return elementAt(index)->uri;
	}

	virtual const XMLCh* getLocalName(const /*unsigned int*/XMLSize_t index) const {
		return elementAt(index)->localPart;
	}

	virtual const XMLCh* getQName(const /*unsigned int*/XMLSize_t index) const {
		return elementAt(index)->qName;
	}

	virtual const XMLCh* getType(const /*unsigned int*/XMLSize_t index) const {
		return elementAt(index)->attrType;
	}

	virtual const XMLCh* getValue(const /*unsigned int*/XMLSize_t index) const {
		return elementAt(index)->value;
	}

	virtual int getIndex(const XMLCh* const uri, const XMLCh* const localPart) const {
		for (unsigned int i=0; i<size(); i++)
			if (XMLString::equals(elementAt(i)->uri, uri) && XMLString::equals(elementAt(i)->localPart, localPart))
				return i;
		return -1;
	}

	virtual int getIndex(const XMLCh* const qName) const {
		for (unsigned int i=0; i<size(); i++)
			if (XMLString::equals(elementAt(i)->qName, qName))
				return i;
		return -1;
	}

	// MT - because new xerces-c but not used
	virtual bool getIndex(const XMLCh* const qName,
                          XMLSize_t& index) const { return false; } ;

	virtual bool getIndex(const XMLCh* const uri,
                          const XMLCh* const localPart,
                          XMLSize_t& index) const { return false; } ;
	// MT

	virtual const XMLCh* getType(const XMLCh* const uri,
			const XMLCh* const localPart) const {
		for (unsigned int i=0; i<size(); i++)
			if (XMLString::equals(elementAt(i)->uri, uri) && XMLString::equals(elementAt(i)->localPart, localPart))
				return elementAt(i)->attrType;
		return NULL;
	}

	virtual const XMLCh* getType(const XMLCh* const qName) const {
		for (unsigned int i=0; i<size(); i++)
			if (XMLString::equals(elementAt(i)->qName, qName))
				return elementAt(i)->attrType;
		return NULL;
	}

	virtual const XMLCh* getValue(const XMLCh* const uri,
			const XMLCh* const localPart) const {
		for (unsigned int i=0; i<size(); i++)
			if (XMLString::equals(elementAt(i)->uri, uri) && XMLString::equals(elementAt(i)->localPart, localPart))
				return elementAt(i)->value;
		return NULL;
	}

	virtual const XMLCh* getValue(const XMLCh* const qName) const {
		for (unsigned int i=0; i<size(); i++)
			if (XMLString::equals(elementAt(i)->qName, qName))
				return elementAt(i)->value;
		return NULL;
	}

};

// This is a simple class that lets us do easy (though not terribly efficient)
// trancoding of XMLCh data to local code page for display.
class StrX {
public:

	// Constructors and Destructor

	StrX(const XMLCh* const toTranscode) {
		// Call the private transcoding method
		fLocalForm = XMLString::transcode(toTranscode);
	}

	~StrX() {
		XMLString::release(&fLocalForm);
	}

	// Getter methods

	const char* localForm() const {
		return fLocalForm;
	}

private:

	// This is the local code page form of the string
	char* fLocalForm;

};

inline XERCES_STD_QUALIFIER ostream& operator<<(XERCES_STD_QUALIFIER ostream& target, const StrX& toDump)
{
	target << toDump.localForm();
	return target;
}

class SAX2SortAttributesFilter : public SAX2XMLFilterImpl {

public:

	SAX2SortAttributesFilter(SAX2XMLReader* parent) :
		SAX2XMLFilterImpl(parent) {
	}

	~SAX2SortAttributesFilter() {
	}

	// Implementations of the SAX2XMLFilter interface	
	void startElement(const XMLCh* const uri, const XMLCh* const localname,
			const XMLCh* const qname, const Attributes& attributes) {
		AttrList sortedList((int)attributes.getLength()); // MT
		for (unsigned int i=0; i<attributes.getLength(); i++) {
			unsigned int j;
			for (j=0; j<sortedList.getLength(); j++) {
				if (XMLString::compareString(sortedList.elementAt(j)->qName, attributes.getQName(i))>=0)
					break;
			}
			Attr* pClone=new Attr;
			pClone->qName = attributes.getQName(i);
			pClone->uri = attributes.getURI(i);
			pClone->localPart = attributes.getLocalName(i);
			pClone->value = attributes.getValue(i);
			pClone->attrType = attributes.getType(i);
			sortedList.insertElementAt(pClone, j);
		}
		SAX2XMLFilterImpl::startElement(uri, localname, qname, sortedList);
	}

};

class Handlers : public DefaultHandler, private XMLFormatTarget {

private:

	// Used by the parser
	XMLFormatter fFormatter;
	bool fExpandNS;
	int state;

	// The parsed XLK file
	std::string xmlFileName;

	// The IFSData object
	IFSData *data;

public:

	Handlers(const char* const encodingName,
			const XMLFormatter::UnRepFlags unRepFlags,
			const bool expandNamespaces, std::string s, IFSData &d) :
		fFormatter(encodingName, 0, this, XMLFormatter::NoEscapes, unRepFlags),
				fExpandNS(expandNamespaces) {
		xmlFileName = s;
		data = &d;
	}

	~Handlers() {
	}

	// Overrides of the output formatter target interface

	void writeChars(const XMLByte* const /* toWrite */) {
	}

	void writeChars(const XMLByte* const toWrite, const /*unsigned int*/XMLSize_t count,
			XMLFormatter* const /* formatter */) {
		XERCES_STD_QUALIFIER cout.write((char *) toWrite, (int) count);
		XERCES_STD_QUALIFIER cout.flush();
	}

	// Overrides of the SAX ErrorHandler interface

	void error(const SAXParseException& e) {
		XERCES_STD_QUALIFIER cerr << "\nError at file " << StrX(e.getSystemId())
		<< ", line " << e.getLineNumber()
		<< ", char " << e.getColumnNumber()
		<< "\n  Message: " << StrX(e.getMessage()) << XERCES_STD_QUALIFIER endl;
	}

	void fatalError(const SAXParseException& e) {
		XERCES_STD_QUALIFIER cerr << "\nFatal Error at file " << StrX(e.getSystemId())
		<< ", line " << e.getLineNumber()
		<< ", char " << e.getColumnNumber()
		<< "\n  Message: " << StrX(e.getMessage()) << XERCES_STD_QUALIFIER endl;
	}

	void warning(const SAXParseException& e) {
		XERCES_STD_QUALIFIER cerr << "\nWarning at file " << StrX(e.getSystemId())
		<< ", line " << e.getLineNumber()
		<< ", char " << e.getColumnNumber()
		<< "\n  Message: " << StrX(e.getMessage()) << XERCES_STD_QUALIFIER endl;
	}

	// Overrides of the SAX DTDHandler interface

	void unparsedEntityDecl(const XMLCh* const /* name */
	, const XMLCh* const /* publicId */
	, const XMLCh* const /* systemId */
	, const XMLCh* const /* notationName */) {
	}

	void notationDecl(const XMLCh* const /* name */
	, const XMLCh* const /* publicId */
	, const XMLCh* const /* systemId */) {
	}

	// Overrides of the SAX DocumentHandler interface

	void characters(const XMLCh* const chars, const unsigned int length) {
	}

	void ignorableWhitespace(const XMLCh* const chars, const unsigned int length) {
		fFormatter.formatBuf(chars, length, XMLFormatter::NoEscapes);
	}

	void processingInstruction(const XMLCh* const target,
			const XMLCh* const data) {
		fFormatter << XMLFormatter::NoEscapes << gStartPI << target;
		if (data)
			fFormatter << chSpace << data;
		fFormatter << XMLFormatter::NoEscapes << gEndPI;
	}

	// General processing methods

	void startDocument();

	void endDocument();

	void startElement(const XMLCh* const uri, const XMLCh* const localname,
			const XMLCh* const qname, const Attributes& attributes);

	void endElement(const XMLCh* const uri, const XMLCh* const localname,
			const XMLCh* const qname);

	void startAttributes(const Attributes& attributes);

	// Process the IndexedFaceSet colorPerVertex attribute
	void parseColorPerVertex(const char* val);

	// Process the IndexedFaceSet coordIndex attribute (faces)
	void parseCoordIndex(const char* val);

	// Process the Coordinate point attribute (vertices)
	void parsePoint(const char* val);

	// Process the COLOR color attribute
	void parseColor(const char* val);

};

// X3DMeshExtractor class

class X3DMeshExtractor {

public:

	X3DMeshExtractor() {
	}

	int load(std::string xmlFileName, IFSData &data);

};

#endif /*X3DMESHEXTRACTOR_H_*/
