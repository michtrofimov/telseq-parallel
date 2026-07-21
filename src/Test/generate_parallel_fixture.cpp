#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"

static BamTools::BamAlignment make_alignment(
    const std::string& name,
    const int refID,
    const int position,
    const std::string& sequence,
    const bool mapped,
    const bool duplicate)
{
    BamTools::BamAlignment alignment;
    alignment.Name = name;
    alignment.RefID = refID;
    alignment.Position = position;
    alignment.QueryBases = sequence;
    alignment.Qualities = std::string(sequence.size(), 'I');
    alignment.MateRefID = -1;
    alignment.MatePosition = -1;
    alignment.InsertSize = 0;
    alignment.MapQuality = mapped ? 60 : 0;
    alignment.SetIsMapped(mapped);
    alignment.SetIsDuplicate(duplicate);

    if (mapped) {
        alignment.CigarData.push_back(
            BamTools::CigarOp(
                'M',
                static_cast<uint32_t>(sequence.size())));
    }

    alignment.AddTag("RG", "Z", std::string("rg1"));
    return alignment;
}

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 4) {
        std::cerr << "Usage: " << argv[0]
                  << " OUTPUT.bam [READS_PER_REFERENCE [READ_LENGTH]]\n";
        return 2;
    }

    const std::string bamPath = argv[1];
    const int readsPerReference = argc >= 3 ? std::atoi(argv[2]) : 20;
    const int readLength = argc >= 4 ? std::atoi(argv[3]) : 100;
    if (readsPerReference < 1 || readLength < 100 || readLength > 100000) {
        std::cerr << "READS_PER_REFERENCE must be positive and READ_LENGTH "
                  << "must be between 100 and 100000\n";
        return 2;
    }

    const int referenceCount = 64;
    const int recordSpacing = readLength + 20;
    const int referenceLength = std::max(
        100000,
        10000 + (readsPerReference + 1) * recordSpacing);
    BamTools::RefVector references;
    std::ostringstream header;
    header << "@HD\tVN:1.0\tSO:coordinate\n";

    for (int refID = 0; refID < referenceCount; ++refID) {
        const std::string name = "contig" + std::to_string(refID);
        references.push_back(
            BamTools::RefData(name, referenceLength));
        header << "@SQ\tSN:" << name
               << "\tLN:" << referenceLength << "\n";
    }
    header << "@RG\tID:rg1\tSM:sample\tLB:library\n";

    BamTools::BamWriter writer;
    if (!writer.Open(bamPath, header.str(), references)) {
        std::cerr << "Could not create " << bamPath << ": "
                  << writer.GetErrorString() << "\n";
        return 1;
    }

    const std::string ordinary(readLength, 'A');
    std::string telomeric =
        "TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    telomeric.resize(readLength, 'A');

    uint64_t physicalRecords = 0;
    uint64_t mappedRecords = 0;
    uint64_t noCoordinateRecords = 0;

    for (int refID = 0; refID < referenceCount; ++refID) {
        // Eight references are deliberately empty. Empty index regions must
        // not be treated as seek failures.
        if (refID % 8 == 7) {
            continue;
        }

        for (int i = 0; i < readsPerReference; ++i) {
            BamTools::BamAlignment alignment = make_alignment(
                "mapped-" + std::to_string(refID) + "-" +
                    std::to_string(i),
                refID,
                100 + i * recordSpacing,
                i == 0 ? telomeric : ordinary,
                true,
                refID == 0 && i == 1);

            if (i % 2 == 1) {
                alignment.SetIsReverseStrand(true);
            }
            if (refID == 1 && i == 2) {
                alignment.SetIsPrimaryAlignment(false);
            }
            if (refID == 2 && i == 3) {
                alignment.SetIsFailedQC(true);
            }

            if (!writer.SaveAlignment(alignment)) {
                std::cerr << writer.GetErrorString() << "\n";
                return 1;
            }
            physicalRecords += 1;
            mappedRecords += 1;
        }

        // This record is unmapped by flag, but still has a coordinate. It must
        // be returned by the indexed reference scan, not by the no-coordinate
        // compatibility branch.
        if (refID == 0) {
            BamTools::BamAlignment positionedUnmapped = make_alignment(
                "positioned-unmapped",
                refID,
                5000,
                ordinary,
                false,
                false);
            if (!writer.SaveAlignment(positionedUnmapped)) {
                std::cerr << writer.GetErrorString() << "\n";
                return 1;
            }
            physicalRecords += 1;
        }
    }

    for (int i = 0; i < 8; ++i) {
        BamTools::BamAlignment alignment = make_alignment(
            "no-coordinate-" + std::to_string(i),
            -1,
            -1,
            i == 7 ? telomeric : ordinary,
            false,
            i == 7);
        if (!writer.SaveAlignment(alignment)) {
            std::cerr << writer.GetErrorString() << "\n";
            return 1;
        }
        physicalRecords += 1;
        noCoordinateRecords += 1;
    }

    writer.Close();

    BamTools::BamReader reader;
    if (!reader.Open(bamPath)) {
        std::cerr << "Could not reopen fixture: "
                  << reader.GetErrorString() << "\n";
        return 1;
    }
    if (!reader.CreateIndex(BamTools::BamIndex::STANDARD)) {
        std::cerr << "Could not index fixture: "
                  << reader.GetErrorString() << "\n";
        return 1;
    }
    reader.Close();

    std::cerr << "Created " << bamPath << " with "
              << referenceCount << " references, "
              << physicalRecords << " physical records, "
              << mappedRecords << " mapped records, and "
              << noCoordinateRecords << " no-coordinate records\n";
    return 0;
}
