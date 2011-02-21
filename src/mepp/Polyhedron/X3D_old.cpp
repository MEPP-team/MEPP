#include "X3D_old.h"

// Write the data on disk
//#define WRITE_ON_DISK


// STATIC DATA //

static const char* encodingName = "LATIN1";
static XMLFormatter::UnRepFlags unRepFlags = XMLFormatter::UnRep_CharRef;
static SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Never;
static bool expandNamespaces= false;
static bool doNamespaces = false;
static bool doSchema = false;
static bool schemaFullChecking = false;
static bool namespacePrefixes = false;
static bool sortAttributes = false;

// We start to parse the document
void Handlers::startDocument() {
	state = NONE;
	data->colorMode = NO_COLOR;

#ifdef WRITE_ON_DISK

	data->openFiles(xmlFileName);

#endif

}

// We have finish the parsing but some data are not exported yet
void Handlers::endDocument() {

#ifdef WRITE_ON_DISK

	data->infoFile << "NumberOfVertices=" << data->vertex.size() << XERCES_STD_QUALIFIER endl;
	data->infoFile << "NumberOfFaces=" << data->face.size() << XERCES_STD_QUALIFIER endl;
	data->infoFile << "NumberOfColors=" << data->color.size() << XERCES_STD_QUALIFIER endl;
	data->closeFiles();

#endif

}

// We enter a new node : change the parser state
void Handlers::startElement(const XMLCh* const uri,
		const XMLCh* const localname, const XMLCh* const qname,
		const Attributes& attributes) {
	if (fExpandNS ) {
		if (XMLString::compareIString(uri, XMLUni::fgZeroLenString) != 0)
			fFormatter << uri << chColon;
		fFormatter << localname ;
	} else {
		char* str = XMLString::transcode(qname);
		if (XMLString::equals(str, "Shape")) {
			state = SHAPE;
		} else if (state == SHAPE&& XMLString::equals(str, "IndexedFaceSet")) {
			state = INDEXEDFACESET;
			startAttributes(attributes);
		} else if (state == INDEXEDFACESET&& XMLString::equals(str,
				"Coordinate")) {
			state = COORDINATE;
			startAttributes(attributes);
		} else if (state == INDEXEDFACESET&& XMLString::equals(str, "Color")) {
			state = COLOR;
			startAttributes(attributes);
		}
		XMLString::release(&str);
	}
}

// We leave a node : change the parser state
void Handlers::endElement(const XMLCh* const uri, const XMLCh* const localname,
		const XMLCh* const qname) {
	if (fExpandNS ) {
		if (XMLString::compareIString(uri, XMLUni::fgZeroLenString) != 0)
			fFormatter << uri << chColon;
		fFormatter << localname;
		fFormatter << localname << chCloseAngle;
	} else {
		char* str = XMLString::transcode(qname);
		if (XMLString::equals(str, "Shape")) {
			state = NONE;
		} else if (XMLString::equals(str, "IndexedFaceSet")) {
			state = SHAPE;
		} else if (state == COORDINATE&& XMLString::equals(str, "Coordinate")) {
			state = INDEXEDFACESET;
		} else if (state == COLOR&& XMLString::equals(str, "Color")) {
			state = INDEXEDFACESET;
		}
		XMLString::release(&str);
	}

}

// We meet a node attribute : process it !
void Handlers::startAttributes(const Attributes& attributes) {

	int len = (int)attributes.getLength(); // MT
	for (int index = 0; index < len; index++) {
		if (fExpandNS ) {
			if (XMLString::compareIString(attributes.getURI(index),
					XMLUni::fgZeroLenString) != 0)
				fFormatter << attributes.getURI(index) << chColon;
			fFormatter << attributes.getLocalName(index) ;
		} else {
			char* str = XMLString::transcode(attributes.getQName(index));
			char* val= NULL;
			if (state == INDEXEDFACESET&& XMLString::equals(str,
					"colorPerVertex")) {
				parseColorPerVertex(val
						= XMLString::transcode(attributes.getValue(index)));
			} else if (state == INDEXEDFACESET&& XMLString::equals(str,
					"coordIndex")) {
				parseCoordIndex(val
						= XMLString::transcode(attributes.getValue(index)));
			} else if (state == COORDINATE&& XMLString::equals(str, "point")) {
				parsePoint(val
						= XMLString::transcode(attributes.getValue(index)));
			} else if (state == COLOR&& XMLString::equals(str, "color")) {
				parseColor(val
						= XMLString::transcode(attributes.getValue(index)));
			}
			XMLString::release(&str);
			if (val != NULL)
				XMLString::release(&val);
		}
	}
}

// process the IndexedFaceSet colorPerVertex attribute
void Handlers::parseColorPerVertex(const char* val) {
	if (data->colorMode == NO_COLOR) {
		if (strcmp(val, "true") == 0) {
			data->colorMode = COLOR_PER_VERTEX;
		} else if (strcmp(val, "false") == 0) {
			data->colorMode = COLOR_PER_FACE;
		}

#ifdef WRITE_ON_DISK

		data->infoFile << "ColorPerVertex=" << data->colorMode << XERCES_STD_QUALIFIER endl;

#endif

	} else {
		XERCES_STD_QUALIFIER cerr << "Different color modes in the same file is *NOT* handled !" << XERCES_STD_QUALIFIER endl;
		XERCES_STD_QUALIFIER cerr << "The parsing job goes on but it may produce invalid export !" << XERCES_STD_QUALIFIER endl;
	}
}

// Process the IndexedFaceSet coordIndex attribute (faces)
void Handlers::parseCoordIndex(const char* val) {
	std::istringstream str(val);

	int index = -2; // we need a value less than -1
	int indexOffset = (int)data->vertex.size(); // MT
	std::vector <int> lastFace;

	while (!str.eof()) {
		str >> index;
		if (str.good()) {
			if (index > -1) {
				lastFace.push_back(indexOffset + index);
			} else if (index == -1) {
				if (lastFace.size() > 0) {

					data->face.push_back(lastFace);
					lastFace.clear();

#ifdef WRITE_ON_DISK

					data->write_last_face();

#endif
				}
			}
		} else {
			str.clear(std::ios::goodbit);
			str.get();
		}
	}

	// There is no delimiter before EOF
	if (lastFace.size() > 0) {
		data->face.push_back(lastFace);
		lastFace.clear();

#ifdef WRITE_ON_DISK

		data->write_last_face();

#endif

	}
}

// Process the Coordinate point attribute (vertices)
void Handlers::parsePoint(const char* val) {
	std::istringstream str(val);
	std::vector <float> coord;
	float tmp;
	int i = 0;

	while (!str.eof()) {
		str >> tmp;

		if (str.good()) {
			coord.push_back(tmp);

			if (i == 2) {
				data->vertex.push_back(coord);
				coord.clear();
				i = 0;

#ifdef WRITE_ON_DISK

				data->write_last_vertex();
#endif
			} else {
				i++;
			}
		} else {
			str.clear(std::ios::goodbit);
			str.get();
		}
	}

	// There is no delimiter before EOF
	if (i == 2) {
		coord.push_back(tmp);
		data->vertex.push_back(coord);
		coord.clear();

#ifdef WRITE_ON_DISK

		data->write_last_vertex();
#endif

	}

}

// Process the COLOR color attribute
void Handlers::parseColor(const char* val) {
	std::istringstream str(val);
	std::vector <float> values;
	float tmp;
	int i = 0;

	while (!str.eof()) {
		str >> tmp;
		if (str.good()) {
			values.push_back(tmp);

			if (i == 2) {
				data->color.push_back(values);
				values.clear();
				i = 0;

#ifdef WRITE_ON_DISK

				data->write_last_color();
#endif
			} else {
				i++;
			}
		} else {
			str.clear(std::ios::goodbit);
			str.get();
		}
	}

	// There is no delimiter before EOF
	if (i == 2) {
		values.push_back(tmp);
		data->color.push_back(values);
		values.clear();

#ifdef WRITE_ON_DISK

		data->write_last_color();
#endif

	}

}
;

int X3DMeshExtractor::load(std::string xmlFileName, IFSData &data) {
	// Initialize the XML4C2 system
	try
	{
		XMLPlatformUtils::Initialize();
	}
	catch (const XMLException& toCatch)
	{
		XERCES_STD_QUALIFIER cerr << "Error during initialization! :\n"
		<< StrX(toCatch.getMessage()) << XERCES_STD_QUALIFIER endl;
		return 1;
	}

	//  Create a SAX parser object. Then, according to what we were told on
	//  the command line, set it to validate or not.
	SAX2XMLReader* parser;
	SAX2XMLReader* filter = NULL;
	SAX2XMLReader* reader = XMLReaderFactory::createXMLReader();

	if (sortAttributes) {
		parser = filter = new SAX2SortAttributesFilter(reader);
	} else
		parser = reader;

	//  Then, according to what we were told on
	//  the command line, set it to validate or not.
	if (valScheme == SAX2XMLReader::Val_Auto) {
		parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
		parser->setFeature(XMLUni::fgXercesDynamic, true);
	}

	if (valScheme == SAX2XMLReader::Val_Never) {
		parser->setFeature(XMLUni::fgSAX2CoreValidation, false);
	}

	if (valScheme == SAX2XMLReader::Val_Always) {
		parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
		parser->setFeature(XMLUni::fgXercesDynamic, false);
	}

	parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, doNamespaces);
	parser->setFeature(XMLUni::fgXercesSchema, doSchema);
	parser->setFeature(XMLUni::fgXercesSchemaFullChecking, schemaFullChecking);
	parser->setFeature(XMLUni::fgSAX2CoreNameSpacePrefixes, namespacePrefixes);

	// Do NOT load online DTD ! 
	parser->setFeature(XMLUni::fgXercesLoadExternalDTD, false);

	//  Create the handler object and install it as the document and error
	//  handler for the parser. Then parse the file and catch any exceptions
	//  that propogate out
	int errorCount = 0;
	int errorCode = 0;

	try
	{
		Handlers handler(encodingName, unRepFlags, expandNamespaces, xmlFileName, data);
		parser->setContentHandler(&handler);
		parser->setErrorHandler(&handler);
		parser->parse(xmlFileName.c_str());
		errorCount = (int)parser->getErrorCount(); // MT
	}
	catch (const OutOfMemoryException&)
	{
		XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
		errorCode = 5;
	}
	catch (const XMLException& toCatch)
	{
		XERCES_STD_QUALIFIER cerr << "\nAn error occurred\n  Error: "
		<< StrX(toCatch.getMessage())
		<< "\n" << XERCES_STD_QUALIFIER endl;
		errorCode = 4;
	}

	if (errorCode) {
		XMLPlatformUtils::Terminate();
		return errorCode;
	}

	//  Delete the parser itself.  Must be done prior to calling Terminate, below.
	delete reader;
	if (filter)
		delete filter;

	// And call the termination method
	XMLPlatformUtils::Terminate();

	if (errorCount > 0)
		return 4;
	else
		return 0;

}
