#include "utils.h"

bool CONSIDER_RX = false;

barcode stringToBarcode(const string& str) {
    vector<bool> res;
    for(uint i(0);i<str.size();i++){
        switch (str[i]){
            case 'A':res.push_back(false);res.push_back(false);break;
            case 'C':res.push_back(false);res.push_back(true);break;
            case 'G':res.push_back(true);res.push_back(false);break;
            default:res.push_back(true);res.push_back(true);break;
        }
    }
  return res;
}

string barcodeToString(barcode b) {
    string str(b.size()/2, 'N');
    uint j = 0;
    for (uint i(0);i<b.size();i+=2) {
        if (b[i]) {
            if(b[i+1]){
                str[j] = 'T';
            } else {
                str[j] = 'G';
            }
        } else {
            if (b[i+1]) {
                str[j] = 'C';
            } else {
                str[j] = 'A';
            }
        }
        j++;
    }
    return str;
}

vector<string> splitString(string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

BamRegion stringToBamRegion(BamReader& reader, string s) {
	BamRegion r;

	vector<string> t = splitString(s, ":");
    if (t.size() != 2) {
		fprintf(stderr, "Error when attempting to convert %s to a BamRegion.\n", s.c_str());
        exit(EXIT_FAILURE);
	}

    vector<string> p = splitString(t[1], "-");
    if (p.size() != 2) {
        fprintf(stderr, "Error when attempting to convert %s to a BamRegion.\n", s.c_str());
        exit(EXIT_FAILURE);
    }

	int leftID = reader.GetReferenceID(t[0]);
	if (leftID == -1) {
		fprintf(stderr, "Cannot find refence with ID %s.\n", t[0].c_str());
		exit(EXIT_FAILURE);
	}

	BamRegion res(leftID, stoi(p[0]), leftID, stoi(p[1]));
	return res;
}

string convertToSam(const BamAlignment& a, RefVector m_references) {
	ostringstream m_out;

    // tab-delimited
    // <QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> [ <TAG>:<VTYPE>:<VALUE> [...] ]

    // write name & alignment flag
    m_out << a.Name << '\t' << a.AlignmentFlag << '\t';

    // write reference name
    if ((a.RefID >= 0) && (a.RefID < (int)m_references.size()))
        m_out << m_references[a.RefID].RefName << '\t';
    else
        m_out << "*\t";

    // write position & map quality
    m_out << a.Position + 1 << '\t' << a.MapQuality << '\t';

    // write CIGAR
    const std::vector<CigarOp>& cigarData = a.CigarData;
    if (cigarData.empty())
        m_out << "*\t";
    else {
        std::vector<CigarOp>::const_iterator cigarIter = cigarData.begin();
        std::vector<CigarOp>::const_iterator cigarEnd = cigarData.end();
        for (; cigarIter != cigarEnd; ++cigarIter) {
            const CigarOp& op = (*cigarIter);
            m_out << op.Length << op.Type;
        }
        m_out << '\t';
    }

    // write mate reference name, mate position, & insert size
    if (a.IsPaired() && (a.MateRefID >= 0) && (a.MateRefID < (int)m_references.size())) {
        if (a.MateRefID == a.RefID)
            m_out << "=\t";
        else
            m_out << m_references[a.MateRefID].RefName << '\t';
        m_out << a.MatePosition + 1 << '\t' << a.InsertSize << '\t';
    } else
        m_out << "*\t0\t0\t";

    // write sequence
    if (a.QueryBases.empty())
        m_out << "*\t";
    else
        m_out << a.QueryBases << '\t';

    // write qualities
    if (a.Qualities.empty() || (a.Qualities.at(0) == (char)0xFF))
        m_out << '*';
    else
        m_out << a.Qualities;

    // write tag data
    const char* tagData = a.TagData.c_str();
    const std::size_t tagDataLength = a.TagData.length();

    std::size_t index = 0;
    while (index < tagDataLength) {

        // write tag name
        std::string tagName = a.TagData.substr(index, 2);
        m_out << '\t' << tagName << ':';
        index += 2;

        // get data type
        char type = a.TagData.at(index);
        ++index;
        switch (type) {
            case (Constants::BAM_TAG_TYPE_ASCII):
                m_out << "A:" << tagData[index];
                ++index;
                break;

            case (Constants::BAM_TAG_TYPE_INT8):
                // force value into integer-type (instead of char value)
                m_out << "i:" << static_cast<int16_t>(tagData[index]);
                ++index;
                break;

            case (Constants::BAM_TAG_TYPE_UINT8):
                // force value into integer-type (instead of char value)
                m_out << "i:" << static_cast<uint16_t>(tagData[index]);
                ++index;
                break;

            case (Constants::BAM_TAG_TYPE_INT16):
                m_out << "i:" << BamTools::UnpackSignedShort(&tagData[index]);
                index += sizeof(int16_t);
                break;

            case (Constants::BAM_TAG_TYPE_UINT16):
                m_out << "i:" << BamTools::UnpackUnsignedShort(&tagData[index]);
                index += sizeof(uint16_t);
                break;

            case (Constants::BAM_TAG_TYPE_INT32):
                m_out << "i:" << BamTools::UnpackSignedInt(&tagData[index]);
                index += sizeof(int32_t);
                break;

            case (Constants::BAM_TAG_TYPE_UINT32):
                m_out << "i:" << BamTools::UnpackUnsignedInt(&tagData[index]);
                index += sizeof(uint32_t);
                break;

            case (Constants::BAM_TAG_TYPE_FLOAT):
                m_out << "f:" << BamTools::UnpackFloat(&tagData[index]);
                index += sizeof(float);
                break;

            case (Constants::BAM_TAG_TYPE_HEX):  // fall-through
            case (Constants::BAM_TAG_TYPE_STRING):
                m_out << type << ':';
                while (tagData[index]) {
                    m_out << tagData[index];
                    ++index;
                }
                ++index;
                break;
        }

        if (tagData[index] == '\0') break;
    }

    return m_out.str();
}