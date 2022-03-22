#include "utils.h"
#include <regex.h>

SequencingTechnology techno = Undefined;

SequencingTechnology determineSequencingTechnology(const string& barcode) {
    SequencingTechnology res;

    regex_t HaplotaggingReg;
    const char* HaplotaggingPattern = "^A[0-9][0-9]C[0-9][0-9]B[0-9][0-9]D[0-9][0-9]$";
    regex_t stLFRReg;
    const char* stLFRPattern = "^[0-9]+_[0-9]+_[0-9]+$";
    regex_t TELLSeqReg;
    const char* TELLSeqPattern = "^[ACGT]+$";

    int rc;
    size_t nmatch = 1;
    regmatch_t pmatch[1];

    if (0 != (rc = regcomp(&HaplotaggingReg, HaplotaggingPattern, REG_EXTENDED))) {
        throw runtime_error("regcomp: An error occured with pattern " + string(HaplotaggingPattern) + ".");
    }

    if (0 != (rc = regcomp(&stLFRReg, stLFRPattern, REG_EXTENDED))) {
        throw runtime_error("regcomp: An error occured with pattern " + string(stLFRPattern) + ".");
    }

    if (0 != (rc = regcomp(&TELLSeqReg, TELLSeqPattern, REG_EXTENDED))) {
        throw runtime_error("regcomp: An error occured with pattern " + string(TELLSeqPattern) + ".");
    }

    if (0 == (rc = regexec(&TELLSeqReg, barcode.c_str(), nmatch, pmatch, 0))) {
        res = TELLSeq;
    } else if (barcode.substr(barcode.length() - 2) == "-1") {
        res = TenX;
    } else if (0 == (rc = regexec(&HaplotaggingReg, barcode.c_str(), nmatch, pmatch, 0))) {
        res = Haplotagging;
    } else if (0 == (rc = regexec(&stLFRReg, barcode.c_str(), nmatch, pmatch, 0))) {
        res = stLFR;
    } else {
        throw runtime_error("determineSequencingTechnology: Unrecognized sequencing technology. Please make sure your barcodes originate from a compatible technology or are reported as nucleotides in the BX:Z tag.");
    }

    return res;
}

string retrieveNucleotidesContent(const string& barcode) {
    // Determine the sequencing technology according to the barcode, and throw an exception if the technology cannot be determined.
    if (techno == Undefined) {
        techno = determineSequencingTechnology(barcode);
    }

    string res;
    switch (techno) {
        case TenX:
            res = barcode.substr(0, barcode.length() - 2);
            break;
        case TELLSeq:
            res = barcode;
            break;
        case Haplotagging:
            res = barcodes_Haplotagging_A[stoi(barcode.substr(1,2)) - 1] + barcodes_Haplotagging_C[stoi(barcode.substr(4,2)) - 1] + barcodes_Haplotagging_B[stoi(barcode.substr(7,2)) - 1] + barcodes_Haplotagging_D[stoi(barcode.substr(10,2)) - 1];
            break;
        case stLFR:
            {
            vector<string> v = splitString(barcode, "_");
            res = barcodes_stLFR[stoi(v[0]) - 1] + barcodes_stLFR[stoi(v[1]) - 1] + barcodes_stLFR[stoi(v[2]) - 1];
            break;
            }
        default:
            throw runtime_error("retrieveNucleotidesContent: Unexpected error. Please make sure your data is valid and attempt running LRez again.");
    }

    return res;
}

bool isValidBarcode(const string& barcode) {
    if (barcode.empty()) {
        return false;
    }

    if (techno == Undefined) {
        techno = determineSequencingTechnology(barcode);
    }

    bool res;
    switch (techno) {
        case TenX:
            res = (barcode.find("N") == string::npos);
            break;
        case TELLSeq:
            res = (barcode.find("N") == string::npos);
            break;
        case Haplotagging:
            res = (barcode.find("00") == string::npos);
            break;
        case stLFR:
            res = (barcode != "0_0_0");
            break;
        default:
            throw runtime_error("isValidBarcode: Unexpected error. Please make sure your data is valid and attempt running LRez again.");
    }

    return res;
}

barcode stringToBarcode(const string& s) {
    string str = retrieveNucleotidesContent(s);
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

vector<string> extractRegions(string contig, int32_t contigSize, unsigned regionSize) {
    vector<string> res;
    for (int32_t i = 0; i < contigSize; i += regionSize) {
        res.push_back(contig + ":" + to_string(i) + "-" + to_string(i + regionSize - 1));
    }
    //res.push_back(contig + ":" + to_string(contigSize - regionSize + 1) + "-" + to_string(contigSize));

    return res;
}

vector<string> extractRegionsList(BamReader& reader, unsigned regionSize) {
    vector<string> regionsList;

     // Get a vector containing reference sequences data
    RefVector rv = reader.GetReferenceData();

    for (RefData d : rv) {
        int id = reader.GetReferenceID(d.RefName);
        BamAlignment al;
        if (id == -1) {
            throw runtime_error("GetReferenceID: Cannot find reference with ID " + d.RefName + ".");
        }   

        // Only process the chromosome if it has alignments
        if (!reader.SetRegion(id, 0, id, d.RefLength - 1)) {
            throw runtime_error("Error while attempting to jump to region " + d.RefName + ":0-" + to_string(d.RefLength).c_str() + ".");
        }
        if (reader.GetNextAlignment(al)) {
            vector<string> w = extractRegions(d.RefName, d.RefLength, regionSize);
            for (string ww : w) {
                regionsList.push_back(ww);
            }
        }
    }

    return regionsList;
}

BamRegion stringToBamRegion(BamReader& reader, string s) {
	BamRegion r;

	vector<string> t = splitString(s, ":");
    if (t.size() != 2) {
        throw runtime_error("stringToBamRegion: Error when attempting to convert " + s + " to a BamRegion.");
	}

    vector<string> p = splitString(t[1], "-");
    if (p.size() != 2) {
        throw runtime_error("stringToBamRegion: Error when attempting to convert " + s + " to a BamRegion.");
    }

	int leftID = reader.GetReferenceID(t[0]);
	if (leftID == -1) {
        throw runtime_error("GetReferenceID: Cannot find reference with ID " + t[0] + ".");
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
                m_out << "i:" << int(static_cast<int8_t>(tagData[index]));
                ++index;
                break;

            case (Constants::BAM_TAG_TYPE_UINT8):
                // force value into integer-type (instead of char value)
                m_out << "i:" << int(static_cast<uint8_t>(tagData[index]));
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
