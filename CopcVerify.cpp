#include <iostream>
#include <fstream>
#include <deque>
#include <cstring>
#include <cmath>

#include <lazperf/Extractor.hpp>
#include <lazperf/header.hpp>
#include <lazperf/readers.hpp>
#include <lazperf/vlr.hpp>

#include "ProgramArgs.hpp"
#include "Comparison.hpp"

#pragma pack(push, 1)
struct Las
{
    int32_t x;
    int32_t y;
    int32_t z;
    uint16_t intensity;
    uint8_t rn : 4;
    uint8_t nr : 4;
    uint8_t class_flags : 4;
    uint8_t channel : 2;
    uint8_t scan_dir : 1;
    uint8_t eof : 1;
    uint8_t classification;
    uint8_t user_data;
    int16_t scan_angle;
    uint16_t point_source_id;
    double gpstime;
    uint16_t red;
    uint16_t green;
    uint16_t blue;
    uint16_t nir;
};

template<typename T>
struct Range
{
    T low { (std::numeric_limits<T>::max)() };
    T high { (std::numeric_limits<T>::lowest)() };
};

struct RangeCalc
{
    Range<double> x;
    Range<double> y;
    Range<double> z;
    Range<double> gpstime;
};

#pragma pack(pop)

struct VoxelKey
{
    int depth;
    int x;
    int y;
    int z;

    operator std::string () const
    {
        return std::to_string(depth) + "-" + std::to_string(x) + "-" + std::to_string(y) + "-" +
            std::to_string(z);
    }

    VoxelKey parent() const
    {
        return VoxelKey {depth - 1, x >> 1, y >> 1, z >> 1};
    }

    bool operator == (const VoxelKey& other) const
    {
        return depth == other.depth && x == other.x && y == other.y && z == other.z;
    }
};

class Entry
{
public:
    VoxelKey key;
    uint64_t offset;
    int32_t byteSize;
    int32_t pointCount;
    static const int Size = 32;

    bool isChildRef() const
    { return pointCount == -1; }

    static Entry create(std::istream& in)
    {
        Entry e;
        std::vector<char> buf(Entry::Size);
        in.read(buf.data(), buf.size());
        lazperf::LeExtractor s(buf.data(), buf.size());
        s >> e.key.depth >> e.key.x >> e.key.y >> e.key.z;
        s >> e.offset >> e.byteSize >> e.pointCount;
        return e;
    }
};
using Entries = std::deque<Entry>;

class Verifier
{
public:
    Verifier(const std::string& filename, bool dump);

    void run();
    void rawVerify();
private:
    void checkEbVlr(const lazperf::copc_info_vlr& vlr);
    uint64_t traverseTree(const lazperf::copc_info_vlr& vlr, int chunkCount);
    std::vector<char> findVlr(const std::string& user_id, int record_id);
    void verifyRanges(const lazperf::copc_info_vlr& vlr);
    void verifyPage(const Entries& page, const Entries& all);
    void readData(Entry entry);
    Entries processPage(Entries page);
    Entries getHierarchyPage(uint64_t offset, uint64_t size);
    int getChunkCount();

    std::string m_filename;
    bool m_dump;
    std::ifstream m_in;
    int m_ebCount = 0;
    int m_pdrf = 0;
    int m_pointLength = 0;
    lazperf::header14 m_header;
    RangeCalc m_ext;
    uint64_t m_fileSize;
    std::map<int, uint64_t> m_totals;
};

int main(int argc, char *argv[])
{
    using namespace pdal;

    std::string filename;
    bool dump;

    ProgramArgs args;
    args.add("filename", "Filename", filename).setPositional();
    args.add("dump,d", "Dump notable information", dump);

    std::vector<std::string> input;
    for (int i = 1; i < argc; ++i)
        input.push_back(argv[i]);

    try
    {
        args.parse(input);
    }
    catch (const arg_error& err)
    {
        std::cerr << err.what() << "\n";
        std::cerr << "Usage: copcverify <filename>\n";
        return -1;
    }
    Verifier v(filename, dump);
    v.run();
}

Verifier::Verifier(const std::string& filename, bool dump) : m_filename(filename), m_dump(dump)
{
    m_in.open(m_filename.data(), std::ios::binary | std::ios::in);
    m_in.seekg(0, std::ios::end);
    m_fileSize = m_in.tellg();
    m_in.seekg(0);
}

void Verifier::run()
{
    m_header = lazperf::header14::create(m_in);
    lazperf::vlr_header vh = lazperf::vlr_header::create(m_in);
    lazperf::copc_info_vlr copcVlr = lazperf::copc_info_vlr::create(m_in);

    if (m_dump)
    {
        std::cout << "Center X Y Z: " << copcVlr.center_x << " " << copcVlr.center_y <<
            " " << copcVlr.center_z << "\n";
        std::cout << "Root node halfsize: " << copcVlr.halfsize << "\n";
        std::cout << "Root node point spacing: " << copcVlr.spacing << "\n";
        std::cout << "GPS time min/max = " << copcVlr.gpstime_minimum << "/" <<
            copcVlr.gpstime_maximum << "\n";
        std::cout << "Total points (header): " << m_header.point_count << "\n";
        std::cout << "Min X Y Z (header): " << m_header.minx << " " << m_header.miny << " " <<
            m_header.minz << "\n";
        std::cout << "Max X Y Z (header): " << m_header.maxx << " " << m_header.maxy << " " <<
            m_header.maxz << "\n";
    }

    if (m_header.version.major != 1 || m_header.version.minor != 4)
        std::cerr << "Invalid COPC file. Found version " <<
            (int)m_header.version.major << "." << (int)m_header.version.minor <<
            " instead of 1.4\n.";
    if (m_header.header_size != lazperf::header14::Size)
        std::cerr << "Invalid COPC file. Found header size of " << m_header.header_size <<
            " instead of " << lazperf::header14::Size << ".\n";
    if ((m_header.point_format_id & 0x80) == 0)
        std::cerr << "Invalid COPC file. Compression bit (high bit) of point format ID "
            "not set.\n";
    int pdrf = m_header.point_format_id & 0x7F;
    if (pdrf < 6 || pdrf > 8)
        std::cerr << "Invalid COPC file. Point format is " << pdrf << ". Should be 6, 7, or 8.";

    if (vh.user_id != "copc")
        std::cerr << "Invalid COPC VLR header. User ID is '" << vh.user_id <<
            "', not 'copc'.\n";
    if (vh.record_id != 1)
        std::cerr << "Invalid COPC VLR header. Record ID is " << (int)vh.record_id << ", not 1.\n";

    for (int i = 0; i < 11; ++i)
        if (copcVlr.reserved[i])
            std::cerr << "Invalid COPC VLR. COPC field reserved[" << i << "] is " <<
                copcVlr.reserved[i] << ", not 0.\n";

    int baseCount = lazperf::baseCount(m_header.point_format_id);
    if (baseCount == 0)
        throw "Bad point record format '" + std::to_string(m_header.point_format_id) + ".";
    m_ebCount = m_header.point_record_length - baseCount;
    m_pdrf = m_header.pointFormat();
    m_pointLength = m_header.point_record_length;

    int chunkCount = getChunkCount();
    /**
    checkEbVlr(copcVlr);
    **/

    uint64_t readPoints = traverseTree(copcVlr, chunkCount);
    if (m_dump)
    {
        std::cout << "Points per level:\n";
        uint64_t sum = 0;
        for (auto& p : m_totals)
        {
            sum += p.second;
            std::cout << "  " << p.first << ": " << p.second << "\n";
        }
        std::cout << "Total of all levels: " << sum << "\n";
    }

    verifyRanges(copcVlr);

//    rawVerify();
}

int Verifier::getChunkCount()
{
    uint64_t chunkoffset;
    uint32_t chunkcount;

    m_in.seekg(m_header.point_offset);
    m_in.read((char *)&chunkoffset, sizeof(chunkoffset));
    m_in.seekg(chunkoffset + 4); // The first 4 bytes are the version, then the chunk count.
    m_in.read((char *)&chunkcount, sizeof(chunkcount));
    if (chunkcount > (std::numeric_limits<int>::max)())
        std::cout << "Chunk count in chunk table exceeds maximum expected.";
    return (int)chunkcount;
}

void Verifier::rawVerify()
{
    std::vector<char> buf(10000);
    m_in.seekg(0);
    m_in.read(buf.data(), 375);
    if (buf[0] != 'L' || buf[1] != 'A' || buf[2] != 'S' || buf[3] != 'F')
        std::cerr << "Invalid LAS header!\n";
}

/**
void Verifier::checkEbVlr(const lazperf::copc_info_vlr& copcVlr)
{
    if (copcVlr.eb_vlr_offset == 0)
        return;
    uint64_t baseOffset = lazperf::header14::Size + lazperf::laz_vlr::Size + copcVlr.size();

    if (copcVlr.eb_vlr_offset < baseOffset)
        std::cerr << "Invalid COPC VLR. Offset to extra bytes VLR, " << copcVlr.eb_vlr_offset <<
            ", is too small. Should be at least " << baseOffset << ".\n"; 
    uint64_t offset = copcVlr.eb_vlr_offset - lazperf::vlr_header::Size;
    m_in.seekg(offset);
    lazperf::vlr_header vh = lazperf::vlr_header::create(m_in);
    if (vh.user_id != "LASF_Spec")
        std::cerr << "Invalid COPC VLR or extra bytes VLR. Found user ID of '" << vh.user_id <<
            ". Expected 'LASF_Spec'.\n";
    iv (vh.record_id != 4)
        std::cerr << "Invalid COPC VLR or extra bytes VLR. Found record ID of '" << vh.record_id <<
            "'. Expected '4'.\n";
    lazperf::eb_vlr ebVlr = lazperf::eb_vlr::create(m_in);



    if (vlr.eb_vlr_size & 192)
        std::cerr 
}
**/

uint64_t Verifier::traverseTree(const lazperf::copc_info_vlr& vlr, int chunkCount)
{
    Entries allEntries;

    Entries page = getHierarchyPage(vlr.root_hier_offset, vlr.root_hier_size);
    allEntries.insert(allEntries.end(), page.begin(), page.end());
    verifyPage(page, allEntries);
    Entries pageEntries = processPage(page);
    while (pageEntries.size())
    {
        Entry e = pageEntries.front();
        pageEntries.pop_front();
        Entries page = getHierarchyPage(e.offset, e.byteSize);
        allEntries.insert(allEntries.end(), page.begin(), page.end());
        verifyPage(page, allEntries);
        Entries entries = processPage(page);
        pageEntries.insert(pageEntries.end(), entries.begin(), entries.end());
    }

    if (allEntries.size() != chunkCount)
        std::cerr << "Number of hierarchy VLR chunk entries doesn't match chunk table count.\n";

    uint64_t totalPoints;
    for (Entry e : allEntries)
        totalPoints += e.pointCount;
    return totalPoints;
}

void Verifier::verifyPage(const Entries& page, const Entries& all)
{
    auto keyExists = [&all](const VoxelKey& k)
    {
        for (const Entry& e : all)
            if (e.key == k)
                return true;
        return false;
    };

    for (const Entry& e : page)
    {
        VoxelKey parent = e.key.parent();
        if (parent.depth < 0)
            continue;
        if (!keyExists(parent))
            std::cerr << "Hierarchy entry " << std::string(e.key) << " has no parent in "
                "existing hierarchy.\n";
    }
}

// A page is just a list of entries.
Entries Verifier::getHierarchyPage(uint64_t offset, uint64_t size)
{
    Entries page;
    int numEntries = size / Entry::Size;
    if (numEntries * Entry::Size != size)
        std::cerr << "Error in hierarchy page size. Not divisible by entry size.";
    m_in.seekg(offset);
    while (numEntries--)
    {
        page.push_back(Entry::create(m_in));
        Entry& e = page.back();
        if (e.offset + e.byteSize > m_fileSize)
            std::cerr << "Invalid entry - offset(" << e.offset << ") + byteSize(" <<
                e.byteSize << ") is greater than file size(" << m_fileSize << ").\n";
    }
    return page;
}

// Process the page and return any child pages.
Entries Verifier::processPage(Entries page)
{
    Entries children;
    for (Entry entry : page)
    {
        if (entry.isChildRef())
            children.push_back(entry);
        else
        {
            m_totals[entry.key.depth] += entry.pointCount;
            readData(entry);
        }
    }
    return children;
}

// Read the data to make sure it decompresses.
void Verifier::readData(Entry entry)
{
    std::vector<char> buf(entry.byteSize);
    m_in.seekg(entry.offset);
    m_in.read(buf.data(), buf.size());

    /**
    lazperf::reader::chunk_decompressor d(m_pdrf, m_ebCount, buf.data());
    Las l;
    while (entry.pointCount--)
    {
        d.decompress(reinterpret_cast<char *>(&l));

        double x = (l.x * m_header.scale.x) + m_header.offset.x;
        m_ext.x.high = (std::max)(x, m_ext.x.high);
        m_ext.x.low = (std::min)(x, m_ext.x.low);

        double y = (l.y * m_header.scale.y) + m_header.offset.y;
        m_ext.y.high = (std::max)(y, m_ext.y.high);
        m_ext.y.low = (std::min)(y, m_ext.y.low);

        double z = (l.z * m_header.scale.z) + m_header.offset.z;
        m_ext.z.high = (std::max)(z, m_ext.z.high);
        m_ext.z.low = (std::min)(z, m_ext.z.low);

        m_ext.gpstime.high = (std::max)(l.gpstime, m_ext.gpstime.high);
        m_ext.gpstime.low = (std::min)(l.gpstime, m_ext.gpstime.low);
    }
    **/
}

void Verifier::verifyRanges(const lazperf::copc_info_vlr& copcVlr)
{
    // Perhaps we should hardcode epsilon
    auto closeEnough = [](double a, double b, float epsilon)
    {
        using namespace std;
        return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
    };

    std::cerr << std::setprecision(17);
    if (!closeEnough(m_ext.x.low, m_header.minx, .0000001))
        std::cerr << "Minimum X value of " << m_ext.x.low << " doesn't match header minimum of " <<
            m_header.minx << ".\n";
    if (!closeEnough(m_ext.x.high, m_header.maxx, .0000001))
        std::cerr << "Maximum X value of " << m_ext.x.high << " doesn't match header maximum of " <<
            m_header.maxx << ".\n";

    if (!closeEnough(m_ext.y.low, m_header.miny, .0000001))
        std::cerr << "Minimum Y value of " << m_ext.y.low << " doesn't match header minimum of " <<
            m_header.miny << ".\n";
    if (!closeEnough(m_ext.y.high, m_header.maxy, .0000001))
        std::cerr << "Maximum Y value of " << m_ext.y.high << " doesn't match header maximum of " <<
            m_header.maxy << ".\n";

    if (!closeEnough(m_ext.z.low, m_header.minz, .0000001))
        std::cerr << "Minimum Z value of " << m_ext.z.low << " doesn't match header minimum of " <<
            m_header.minz << ".\n";
    if (!closeEnough(m_ext.z.high, m_header.maxz, .0000001))
        std::cerr << "Maximum Z value of " << m_ext.z.high << " doesn't match header maximum of " <<
            m_header.maxz << ".\n";

    if (!closeEnough(m_ext.gpstime.low, copcVlr.gpstime_minimum, .0000001))
        std::cerr << "Minimum GPS time value of " << m_ext.gpstime.low <<
            " doesn't match COPC VLR minimum of " << copcVlr.gpstime_minimum << ".\n";
    if (!closeEnough(m_ext.gpstime.high, copcVlr.gpstime_maximum, .0000001))
        std::cerr << "Maximum GPS time value of " << m_ext.gpstime.high <<
            " doesn't match COPC VLR maximum of " << copcVlr.gpstime_maximum << ".\n";
    std::cerr << std::setprecision(6);
}

std::vector<char> Verifier::findVlr(const std::string& user_id, int record_id)
{
    std::vector<char> buf;

    m_in.seekg(m_header.header_size);
    int count = 0;
    while (count < m_header.vlr_count && m_in.good() && !m_in.eof())    
    {
        lazperf::vlr_header h = lazperf::vlr_header::create(m_in);
        if (h.user_id == user_id && h.record_id == record_id)
        {
            buf.resize(h.data_length);
            m_in.read(buf.data(), h.data_length);
            return buf;
        }
        m_in.seekg(h.data_length, std::ios::cur);
        count++;
    }
    if (m_header.evlr_count == 0 || m_header.evlr_offset == 0)
        return buf;

    m_in.seekg(m_header.evlr_offset);
    count = 0;
    while (count < m_header.evlr_count && m_in.good() && !m_in.eof())    
    {
        lazperf::evlr_header h = lazperf::evlr_header::create(m_in);
        if (h.user_id == user_id && h.record_id == record_id)
        {
            buf.resize(h.data_length);
            m_in.read(buf.data(), h.data_length);
            return buf;
        }
        m_in.seekg(h.data_length, std::ios::cur);
        count++;
    }
    return buf;
}