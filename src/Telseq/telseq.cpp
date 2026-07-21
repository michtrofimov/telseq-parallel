
#include <iostream>
#include <istream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <vector>
#include <istream>
#include <sstream>
#include <set>
#include <map>
#include <atomic>
#include <exception>
#include <memory>
#include <mutex>
#include <thread>
#include <unistd.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "telseq.h"
#include "Timer.h"
#include "math.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "Util.h"
#include "prettyprint.h"

//
// Getopt
//

#define PROGRAM_BIN "telseq"
#define AUTHOR "Zhihao Ding"

static const char *TELSEQ_VERSION_MESSAGE =
"TelSeq Version " PACKAGE_VERSION "\n"
"Written by " AUTHOR ".\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n" ;

static const char *TELSEQ_USAGE_MESSAGE =
"Program: " PACKAGE_NAME "\n"
"Version: " PACKAGE_VERSION "\n"
"Contact: " AUTHOR " [" PACKAGE_BUGREPORT "]\n\n"
"Usage: " PROGRAM_BIN " [OPTION] <in.1.bam> <in.2.bam> <...> \n"
"Scan BAM and estimate telomere length. \n"
"   <in.bam>                 one or more BAM files to be analysed. File names can also be passed from a pipe, \n "
"                            with each row containing one BAM path.\n"
"   -f, --bamlist=STR        a file that contains a list of file paths of BAMs. It should has only one column, \n"
"                            with each row a BAM file path. -f has higher priority than <in.bam>. When specified, \n"
"                            <in.bam> are ignored.\n"
"   -o, --output_dir=STR     output file for results. Ignored when input is from stdin, in which case output will be stdout. \n"
"   -H                       remove header line, which is printed by default.\n"
"   -h                       print the header line only. The text can be used to attach to result files, useful\n"
"                            when the headers of the result files are suppressed. \n"
"   -m                       merge read groups by taking a weighted average across read groups of a sample, weighted by \n"
"                            the total number of reads in read group. Default is to output each readgroup separately.\n"
"   -u                       ignore read groups. Treat all reads in BAM as if they were from a same read group.\n"
"   -t, --threads=INT        number of threads for one coordinate-sorted, indexed BAM. default = 1.\n"
"                            Values greater than 1 reserve one compatibility scanner; remaining workers use the index.\n"
"   -k                       threshold of the amount of TTAGGG/CCCTAA repeats in read for a read to be considered telomeric. default = 7.\n"
"\nTesting functions\n------------\n"
"   -r                       read length. default = 100\n"
"   -z                       use user specified pattern for searching [ATGC]*.\n"
"   -e, --exomebed=STR       specifiy exome regions in BED format. These regions will be excluded \n"
"   -w,                      consider BAMs in the speicfied bamlist as one single BAM. This is useful when \n"
"                            the initial alignemt is separated for some reason, such as one for mapped and one for ummapped reads. \n"
"   --help                   display this help and exit\n"

"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static StringVector bamlist;
    static std::string outputfile = "";
    static std::string exomebedfile = "";
    static std::map< std::string, std::vector<range> > exomebed;
    static bool writerheader = true;
    static bool mergerg = false;
    static bool ignorerg = false;
    static bool onebam = false; // whether to consider all bams as one bam
    static unsigned int threads = 1;
    static int tel_k= ScanParameters::TEL_MOTIF_CUTOFF;
    static std::string unknown = "UNKNOWN";
    static std::string PATTERN;
    static std::string PATTERN_REV;

}

static const char* shortopts = "f:o:k:z:e:r:p:t:Hhvmuw";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bamlist",		optional_argument, NULL, 'f' },
    { "output-dir",		optional_argument, NULL, 'o' },
    { "exomebed",		optional_argument, NULL, 'e' },
    { "threads",            required_argument, NULL, 't' },
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

// combine counts in two ScanResults objects
void add_results(ScanResults& x, ScanResults& y){

    x.numTotal += y.numTotal;
    x.numMapped += y.numMapped;
    x.numDuplicates += y.numMapped;
    x.n_exreadsExcluded += y.n_exreadsExcluded;
    x.n_exreadsChrUnmatched += y.n_exreadsChrUnmatched;
    x.n_totalunfiltered += y.n_totalunfiltered;

    for (std::size_t j = 0, max = x.telcounts.size(); j != max; ++j){
        x.telcounts[j] +=y.telcounts[j];
    }
    for (std::size_t k = 0, max = x.gccounts.size(); k != max; ++k){
        x.gccounts[k] += y.gccounts[k];
    }

}

typedef std::map<std::string, ScanResults> ResultMap;
typedef std::map<std::string, std::vector<range>::iterator> ExomeSearchHints;

struct ReferenceTask
{
    int refID;
    int refLength;
};

static std::mutex parallelLogMutex;

struct SamFileDeleter
{
    void operator()(samFile* file) const
    {
        if (file != NULL) {
            sam_close(file);
        }
    }
};

struct SamHeaderDeleter
{
    void operator()(bam_hdr_t* header) const
    {
        if (header != NULL) {
            bam_hdr_destroy(header);
        }
    }
};

struct HtsIndexDeleter
{
    void operator()(hts_idx_t* index) const
    {
        if (index != NULL) {
            hts_idx_destroy(index);
        }
    }
};

struct HtsIteratorDeleter
{
    void operator()(hts_itr_t* iterator) const
    {
        if (iterator != NULL) {
            hts_itr_destroy(iterator);
        }
    }
};

struct BamRecordDeleter
{
    void operator()(bam1_t* record) const
    {
        if (record != NULL) {
            bam_destroy1(record);
        }
    }
};

typedef std::unique_ptr<samFile, SamFileDeleter> SamFilePtr;
typedef std::unique_ptr<bam_hdr_t, SamHeaderDeleter> SamHeaderPtr;
typedef std::unique_ptr<hts_idx_t, HtsIndexDeleter> HtsIndexPtr;
typedef std::unique_ptr<hts_itr_t, HtsIteratorDeleter> HtsIteratorPtr;
typedef std::unique_ptr<bam1_t, BamRecordDeleter> BamRecordPtr;

// This reducer is intentionally separate from add_results().  The latter is
// part of the legacy -w path and its behaviour is kept unchanged for output
// compatibility with the original TelSeq release.
static void add_thread_results(ScanResults& destination, const ScanResults& source)
{
    destination.numTotal += source.numTotal;
    destination.numMapped += source.numMapped;
    destination.numDuplicates += source.numDuplicates;
    destination.n_exreadsExcluded += source.n_exreadsExcluded;
    destination.n_exreadsChrUnmatched += source.n_exreadsChrUnmatched;
    destination.n_totalunfiltered += source.n_totalunfiltered;

    for (std::size_t i = 0; i < destination.telcounts.size(); ++i) {
        destination.telcounts[i] += source.telcounts[i];
    }
    for (std::size_t i = 0; i < destination.gccounts.size(); ++i) {
        destination.gccounts[i] += source.gccounts[i];
    }
}

static void print_parallel_warning(const std::string& message)
{
    std::lock_guard<std::mutex> lock(parallelLogMutex);
    std::cerr << message;
}

static bool convert_hts_alignment(
    const bam1_t* source,
    BamTools::BamAlignment& destination,
    std::string& error)
{
    if (source == NULL) {
        error = "HTSlib returned a null BAM alignment";
        return false;
    }

    destination = BamTools::BamAlignment();
    destination.RefID = source->core.tid;
    destination.Position = source->core.pos;
    destination.Length = source->core.l_qseq;
    destination.AlignmentFlag = source->core.flag;

    static const char sequenceLookup[] = "=ACMGRSVTWYHKDBN";
    const uint8_t* encodedSequence = bam_get_seq(source);
    destination.QueryBases.resize(source->core.l_qseq);
    for (int i = 0; i < source->core.l_qseq; ++i) {
        destination.QueryBases[i] =
            sequenceLookup[bam_seqi(encodedSequence, i)];
    }

    const uint8_t* readGroupTag = bam_aux_get(source, "RG");
    if (readGroupTag != NULL) {
        const char* readGroup = bam_aux2Z(readGroupTag);
        if (readGroup == NULL) {
            error = "HTSlib found a non-string RG tag";
            return false;
        }
        if (!destination.AddTag("RG", "Z", std::string(readGroup))) {
            error = "could not copy RG tag from HTSlib alignment";
            return false;
        }
    }

    return true;
}

// Scan one alignment into a worker-local result map. This mirrors the body of
// the original serial scan loop, but never touches shared counters.
static void scan_alignment_parallel(
    BamTools::BamAlignment& record,
    const bool rggroups,
    const bool isExome,
    ResultMap& resultmap,
    ExomeSearchHints& lastfound,
    uint64_t& nprocessed)
{
    std::string tag = opt::unknown;
    if (rggroups) {
        if (record.HasTag("RG")) {
            record.GetTag("RG", tag);
        } else {
            std::ostringstream message;
            message << "can't find RG tag for read at position {"
                    << record.RefID << ":" << record.Position << "}\n"
                    << "skip this read\n";
            print_parallel_warning(message.str());
            return;
        }
    }

    ResultMap::iterator result = resultmap.find(tag);
    if (result == resultmap.end()) {
        std::ostringstream message;
        message << "RG tag {" << tag << "} for read at position {"
                << record.RefID << ":" << record.Position
                << "} doesn't exist in BAM header.";
        print_parallel_warning(message.str());
        return;
    }

    if (isExome) {
        range alignmentRange;
        alignmentRange.first = record.Position;
        alignmentRange.second = record.Position + record.Length;
        const std::string chrm = refID2Name(record.RefID);

        if (chrm != "-1") {
            std::map<std::string, std::vector<range> >::iterator chromosome =
                opt::exomebed.find(chrm);
            if (chromosome == opt::exomebed.end()) {
                result->second.n_exreadsChrUnmatched += 1;
            } else {
                const std::vector<range>::iterator end = chromosome->second.end();
                ExomeSearchHints::iterator hint = lastfound.find(chrm);
                if (hint == lastfound.end()) {
                    lastfound[chrm] = chromosome->second.begin();
                }
                std::vector<range>::iterator found =
                    searchRange(lastfound[chrm], end, alignmentRange);
                if (found != end) {
                    result->second.n_exreadsExcluded += 1;
                    lastfound[chrm] = found;
                    return;
                }
            }
        }
    }

    result->second.numTotal += 1;

    if (record.IsMapped()) {
        result->second.numMapped += 1;
    }
    if (record.IsDuplicate()) {
        result->second.numDuplicates += 1;
    }

    const double gc = calcGC(record.QueryBases);
    const int ptn_count =
        countMotif(record.QueryBases, opt::PATTERN, opt::PATTERN_REV);

    if (ptn_count > ScanParameters::TEL_MOTIF_N - 1) {
        return;
    }
    result->second.telcounts[ptn_count] += 1;

    if (gc >= ScanParameters::GC_LOWERBOUND &&
        gc <= ScanParameters::GC_UPPERBOUND) {
        const int idx = floor(
            (gc - ScanParameters::GC_LOWERBOUND) / ScanParameters::GC_BINSIZE);
        assert(idx >= 0 && idx <= ScanParameters::GC_BIN_N - 1);
        if (idx > ScanParameters::GC_BIN_N - 1) {
            std::ostringstream message;
            message << nprocessed << " GC:{" << gc << "} telcounts:{"
                    << ptn_count << "} GC bin index out of bound:" << idx
                    << "\n";
            print_parallel_warning(message.str());
            return;
        }
        result->second.gccounts[idx] += 1;
    }

    nprocessed += 1;
}

static bool scan_htslib_alignment(
    const bam1_t* source,
    const bool rggroups,
    const bool isExome,
    ResultMap& resultmap,
    ExomeSearchHints& lastfound,
    uint64_t& nprocessed,
    std::string& error)
{
    BamTools::BamAlignment converted;
    if (!convert_hts_alignment(source, converted, error)) {
        return false;
    }
    scan_alignment_parallel(
        converted,
        rggroups,
        isExome,
        resultmap,
        lastfound,
        nprocessed);
    return true;
}

static bool header_is_coordinate_sorted(const std::string& headerText)
{
    std::istringstream headerStream(headerText);
    std::string line;
    while (std::getline(headerStream, line)) {
        if (line.compare(0, 3, "@HD") == 0) {
            return line.find("\tSO:coordinate") != std::string::npos;
        }
    }
    return false;
}

static bool scan_bam_parallel(
    const std::string& bamPath,
    BamTools::BamReader& setupReader,
    const ResultMap& emptyResultMap,
    const bool rggroups,
    const bool isExome,
    ResultMap& combinedResults,
    uint64_t& totalScanned)
{
    if (!header_is_coordinate_sorted(setupReader.GetHeaderText())) {
        std::cerr << "Error: -t > 1 requires a BAM whose header declares "
                  << "SO:coordinate.\n";
        return false;
    }

    if (!setupReader.LocateIndex()) {
        std::cerr << "Error: -t > 1 requires a readable BAM index next to "
                  << bamPath << ".\n";
        std::cerr << "BamTools: " << setupReader.GetErrorString() << "\n";
        return false;
    }

    const BamTools::RefVector references = setupReader.GetReferenceData();
    if (references.empty()) {
        std::cerr << "Error: indexed parallel scan found no references in "
                  << bamPath << ".\n";
        return false;
    }

    std::vector<ReferenceTask> tasks;
    tasks.reserve(references.size());
    for (int refID = 0;
         refID < static_cast<int>(references.size());
         ++refID) {
        ReferenceTask task;
        task.refID = refID;
        task.refLength = references[refID].RefLength;
        tasks.push_back(task);
    }

    // Long references start first. Workers fetch another reference as soon as
    // they finish, which balances the small/decoy contigs without splitting a
    // chromosome at overlapping region boundaries.
    std::sort(
        tasks.begin(),
        tasks.end(),
        [](const ReferenceTask& left, const ReferenceTask& right) {
            return left.refLength > right.refLength;
        });

    // One thread uses HTSlib's no-coordinate index iterator for unplaced
    // records and the legacy final-record contribution. The remaining threads
    // process mapped reference tasks through the index.
    const std::size_t mappedWorkerCount =
        std::min<std::size_t>(opt::threads - 1, tasks.size());
    std::vector<ResultMap> workerResults(mappedWorkerCount, emptyResultMap);
    std::vector<uint64_t> workerScanned(mappedWorkerCount, 0);
    std::vector<uint64_t> workerProcessed(mappedWorkerCount, 0);
    std::vector<std::string> workerErrors(mappedWorkerCount);
    ResultMap tailResults = emptyResultMap;
    uint64_t tailScanned = 0;
    uint64_t tailProcessed = 0;
    uint64_t tailIndexedFetched = 0;
    uint64_t tailNoCoordinateFetched = 0;
    uint64_t tailFallbackFetched = 0;
    std::string tailError;
    std::atomic<std::size_t> nextTask(0);
    std::atomic<bool> failed(false);
    std::vector<std::thread> workers;
    workers.reserve(mappedWorkerCount);

    std::cerr << "Indexed parallel scan using " << mappedWorkerCount
              << " mapped-reference workers and 1 HTSlib compatibility scanner "
              << "across " << tasks.size() << " reference tasks\n";

    std::thread tailWorker([&]() {
        try {
            SamFilePtr reader(sam_open(bamPath.c_str(), "rb"));
            if (!reader) {
                tailError = "HTSlib could not open BAM for compatibility scan";
                failed.store(true);
                return;
            }

            SamHeaderPtr header(sam_hdr_read(reader.get()));
            if (!header) {
                tailError = "HTSlib could not read BAM header";
                failed.store(true);
                return;
            }

            HtsIndexPtr index(sam_index_load(reader.get(), bamPath.c_str()));
            if (!index) {
                tailError =
                    "HTSlib could not load a BAI or CSI index for " + bamPath;
                failed.store(true);
                return;
            }

            BamRecordPtr record(bam_init1());
            BamRecordPtr finalRecord(bam_init1());
            if (!record || !finalRecord) {
                tailError = "HTSlib could not allocate BAM records";
                failed.store(true);
                return;
            }

            ExomeSearchHints lastfound;
            bool haveFinalRecord = false;

            HtsIteratorPtr noCoordinateIterator(
                sam_itr_queryi(index.get(), HTS_IDX_NOCOOR, 0, 0));
            if (!noCoordinateIterator) {
                tailError =
                    "HTSlib index does not support direct no-coordinate access";
                failed.store(true);
                return;
            }

            int iteratorStatus = 0;
            while (!failed.load() &&
                   (iteratorStatus = sam_itr_next(
                        reader.get(),
                        noCoordinateIterator.get(),
                        record.get())) >= 0) {
                tailIndexedFetched += 1;
                tailNoCoordinateFetched += 1;
                tailScanned += 1;

                if (!scan_htslib_alignment(
                        record.get(),
                        rggroups,
                        isExome,
                        tailResults,
                        lastfound,
                        tailProcessed,
                        tailError)) {
                    failed.store(true);
                    return;
                }

                if (bam_copy1(finalRecord.get(), record.get()) == NULL) {
                    tailError = "HTSlib could not retain final BAM record";
                    failed.store(true);
                    return;
                }
                haveFinalRecord = true;
            }

            if (iteratorStatus < -1) {
                tailError =
                    "HTSlib failed while reading no-coordinate BAM records";
                failed.store(true);
                return;
            }

            // A coordinate-sorted BAM normally ends with its no-coordinate
            // records, so the last one above is also the final physical BAM
            // record. If that tail is empty, find the highest populated
            // reference and retain its last record. This reads at most one
            // reference rather than performing a full sequential BAM pass.
            if (!haveFinalRecord) {
                for (int refID = header->n_targets - 1;
                     refID >= 0 && !haveFinalRecord;
                     --refID) {
                    uint64_t mapped = 0;
                    uint64_t unmapped = 0;
                    if (hts_idx_get_stat(
                            index.get(), refID, &mapped, &unmapped) == 0 &&
                        mapped + unmapped == 0) {
                        continue;
                    }

                    HtsIteratorPtr referenceIterator(
                        sam_itr_queryi(
                            index.get(),
                            refID,
                            0,
                            header->target_len[refID]));
                    if (!referenceIterator) {
                        tailError =
                            "HTSlib could not query final populated reference";
                        failed.store(true);
                        return;
                    }

                    iteratorStatus = 0;
                    bool referenceHadRecords = false;
                    while (!failed.load() &&
                           (iteratorStatus = sam_itr_next(
                                reader.get(),
                                referenceIterator.get(),
                                record.get())) >= 0) {
                        tailIndexedFetched += 1;
                        tailFallbackFetched += 1;
                        referenceHadRecords = true;
                        if (bam_copy1(finalRecord.get(), record.get()) == NULL) {
                            tailError =
                                "HTSlib could not retain final BAM record";
                            failed.store(true);
                            return;
                        }
                    }

                    if (iteratorStatus < -1) {
                        tailError =
                            "HTSlib failed while locating final BAM record";
                        failed.store(true);
                        return;
                    }
                    haveFinalRecord = referenceHadRecords;
                }
            }

            // The original loop scans its alignment object once after EOF. If
            // the BAM is empty it scans a default object; otherwise it scans
            // the final physical record a second time.
            tailScanned += 1;
            if (haveFinalRecord) {
                if (!scan_htslib_alignment(
                        finalRecord.get(),
                        rggroups,
                        isExome,
                        tailResults,
                        lastfound,
                        tailProcessed,
                        tailError)) {
                    failed.store(true);
                    return;
                }
            } else {
                BamTools::BamAlignment emptyRecord;
                scan_alignment_parallel(
                    emptyRecord,
                    rggroups,
                    isExome,
                    tailResults,
                    lastfound,
                    tailProcessed);
            }
        } catch (const std::exception& error) {
            tailError = error.what();
            failed.store(true);
        } catch (...) {
            tailError = "unknown no-coordinate scanner exception";
            failed.store(true);
        }
    });

    for (std::size_t workerID = 0;
         workerID < mappedWorkerCount;
         ++workerID) {
        workers.push_back(std::thread([&, workerID]() {
            try {
                BamTools::BamReader reader;
                if (!reader.Open(bamPath)) {
                    workerErrors[workerID] =
                        "could not open BAM: " + reader.GetErrorString();
                    failed.store(true);
                    return;
                }
                if (!reader.LocateIndex()) {
                    workerErrors[workerID] =
                        "could not locate BAM index: " + reader.GetErrorString();
                    failed.store(true);
                    reader.Close();
                    return;
                }

                ExomeSearchHints lastfound;
                while (!failed.load()) {
                    const std::size_t taskIndex = nextTask.fetch_add(1);
                    if (taskIndex >= tasks.size()) {
                        break;
                    }

                    const ReferenceTask& task = tasks[taskIndex];
                    const bool positioned = reader.SetRegion(
                        task.refID,
                        0,
                        task.refID,
                        task.refLength);

                    if (!positioned) {
                        workerErrors[workerID] =
                            "could not seek to reference " +
                            references[task.refID].RefName + ": " +
                            reader.GetErrorString();
                        failed.store(true);
                        break;
                    }

                    BamTools::BamAlignment record;
                    while (reader.GetNextAlignment(record)) {
                        workerScanned[workerID] += 1;
                        scan_alignment_parallel(
                            record,
                            rggroups,
                            isExome,
                            workerResults[workerID],
                            lastfound,
                            workerProcessed[workerID]);
                    }
                }
                reader.Close();
            } catch (const std::exception& error) {
                workerErrors[workerID] = error.what();
                failed.store(true);
            } catch (...) {
                workerErrors[workerID] = "unknown worker exception";
                failed.store(true);
            }
        }));
    }

    for (std::size_t i = 0; i < workers.size(); ++i) {
        workers[i].join();
    }
    tailWorker.join();

    if (failed.load()) {
        if (!tailError.empty()) {
            std::cerr << "Error in no-coordinate scanner: "
                      << tailError << "\n";
        }
        for (std::size_t i = 0; i < workerErrors.size(); ++i) {
            if (!workerErrors[i].empty()) {
                std::cerr << "Error in worker " << i + 1 << ": "
                          << workerErrors[i] << "\n";
            }
        }
        return false;
    }

    combinedResults = emptyResultMap;
    totalScanned = tailScanned;
    uint64_t totalProcessed = tailProcessed;
    for (ResultMap::const_iterator source = tailResults.begin();
         source != tailResults.end();
         ++source) {
        ResultMap::iterator destination = combinedResults.find(source->first);
        if (destination == combinedResults.end()) {
            std::cerr << "Error: no-coordinate scanner returned unknown "
                      << "read group " << source->first << "\n";
            return false;
        }
        add_thread_results(destination->second, source->second);
    }

    for (std::size_t workerID = 0;
         workerID < mappedWorkerCount;
         ++workerID) {
        totalScanned += workerScanned[workerID];
        totalProcessed += workerProcessed[workerID];
        for (ResultMap::const_iterator source =
                 workerResults[workerID].begin();
             source != workerResults[workerID].end();
             ++source) {
            ResultMap::iterator destination =
                combinedResults.find(source->first);
            if (destination == combinedResults.end()) {
                std::cerr << "Error: worker returned unknown read group "
                          << source->first << "\n";
                return false;
            }
            add_thread_results(destination->second, source->second);
        }
    }

    std::cerr << "[scan] HTSlib compatibility scanner fetched "
              << tailIndexedFetched << " indexed BAM records: "
              << tailNoCoordinateFetched << " no-coordinate, "
              << tailFallbackFetched << " final-reference fallback; "
              << "full sequential scans: 0\n";
    std::cerr << "[scan] parallel workers processed " << totalProcessed
              << " reads after filters\n";
    return true;
}

// merge results in result list into one
void merge_results_by_readgroup(

    std::vector< std::map<std::string, ScanResults> >& resultlist){
    std::vector< std::map<std::string, ScanResults> > mergedresultslist;
    std::map<std::string, ScanResults> mergedresults;

    for(size_t i =0; i< resultlist.size(); i++){
        auto rmap = resultlist[i];
        for(std::map<std::string, ScanResults>::iterator it= rmap.begin();
                it != rmap.end(); ++it){

            std::string rg = it ->first;
            ScanResults result = it -> second;

            if(mergedresults.find(rg) == mergedresults.end()){
                mergedresults[rg] = result;
            }else{
                add_results(mergedresults[rg], result);
            }
        }
    }
    mergedresultslist.push_back(mergedresults);
    resultlist = mergedresultslist;

}


int scanBam()
{

  std::vector< std::map<std::string, ScanResults> > resultlist;
  bool isExome = opt::exomebedfile.size()==0? false: true;

  std::cout << opt::bamlist << "\n";
  std::cout << opt::bamlist.size() << " BAMs" <<  std::endl;

  for(std::size_t i=0; i<opt::bamlist.size(); i++) {

  // storing results for each read group (RG tag). use
  // read group ID as key.
      std::map<std::string, ScanResults> resultmap;
      // store where the overlap was last found in the case of exome seq
      std::map<std::string, std::vector<range>::iterator> lastfound;
      std::vector<range>::iterator searchhint;

      std::cerr << "Start analysing BAM " << opt::bamlist[i] << "\n";

      // Open the bam files for reading/writing
      BamTools::BamReader* pBamReader = new BamTools::BamReader;

      pBamReader->Open(opt::bamlist[i]);

      // get bam headers
      const BamTools::SamHeader header = pBamReader ->GetHeader();
      bool rggroups=false;

      if(opt::ignorerg){ // ignore read groups
	      std::cerr << "Treat all reads in BAM as if they were from a same sample" << std::endl;
	      ScanResults results;
	      results.sample = opt::unknown;
	      resultmap[opt::unknown]=results;
      }else{
	std::map <std::string, std::string> readgroups;
	std::map <std::string, std::string> readlibs;

	rggroups = header.HasReadGroups();

	if(rggroups){
	  for(BamTools::SamReadGroupConstIterator it = header.ReadGroups.Begin(); it != header.ReadGroups.End();++it){
	    readgroups[it->ID]= it->Sample;
	    if(it->HasLibrary()){
	      readlibs[it->ID] = it -> Library;
	    }else{
	      readlibs[it->ID] = opt::unknown;
	    }
	  }
	  std::cerr<<"Specified BAM has "<< readgroups.size()<< " read groups" << std::endl;

	  for(std::map<std::string, std::string>::iterator it = readgroups.begin(); it != readgroups.end(); ++it){
		  ScanResults results;
		  std::string rgid = it -> first;
		  results.sample = it -> second;
		  results.lib = readlibs[rgid];
		  resultmap[rgid]=results; //results are identified by RG tag.
	  }

	}else{
	  std::cerr << "Warning: can't find RG tag in the BAM header" << std::endl;
	  std::cerr << "Warning: treat all reads in BAM as if they were from a same sample" << std::endl;
	  ScanResults results;
	  results.sample = opt::unknown;
	  results.lib = opt::unknown;
	  resultmap[opt::unknown]=results;
	}
      }

      if(opt::threads > 1){
        ResultMap parallelResults;
        uint64_t parallelTotal = 0;
        const bool parallelOk = scan_bam_parallel(
            opt::bamlist[i],
            *pBamReader,
            resultmap,
            rggroups,
            isExome,
            parallelResults,
            parallelTotal);

        pBamReader->Close();
        delete pBamReader;

        if(!parallelOk){
          return 1;
        }

        resultlist.push_back(parallelResults);
        std::cerr << "[scan] total reads in BAM scanned "
                  << parallelTotal << std::endl;
        std::cerr << "Completed scanning BAM\n";
        continue;
      }

      BamTools::BamAlignment record1;
      bool done = false;

      int nprocessed=0; // number of reads analyzed
      int ntotal=0; // number of reads scanned in bam (we skip some reads, see below)
      while(!done) {
	ntotal ++;
	done = !pBamReader -> GetNextAlignment(record1);
	std::string tag = opt::unknown;
	if(rggroups){
	  // skip reads that do not have read group tag
	  if(record1.HasTag("RG")){
	    record1.GetTag("RG", tag);
	  }else{
	    std::cerr << "can't find RG tag for read at position {" << record1.RefID << ":" << record1.Position << "}" << std::endl;
	    std::cerr << "skip this read" << std::endl;
	    continue;
	  }
	}

	// skip reads with readgroup not defined in BAM header
	if(resultmap.find(tag) == resultmap.end()){
	  std::cerr << "RG tag {" << tag << "} for read at position ";
	  std::cerr << "{" << record1.RefID << ":" << record1.Position << "} doesn't exist in BAM header.";
	  continue;
	}

	// for exome, exclude reads mapped to the exome regions.
	if(isExome){
	  range rg;
	  rg.first = record1.Position;
	  rg.second = record1.Position + record1.Length;
	  std::string chrm =  refID2Name(record1.RefID);

	  if(chrm != "-1"){ // check if overlap exome when the read is mapped to chr1-22, X, Y
	    std::map<std::string, std::vector<range> >::iterator chrmit = opt::exomebed.find(chrm);
	    if(chrmit == opt::exomebed.end()) {
		// unmapped reads can have chr names as a star (*). We also don't consider MT reads.
		resultmap[tag].n_exreadsChrUnmatched +=1;
	    } else {
	      std::vector<range>::iterator itend = opt::exomebed[chrm].end();
	      std::map<std::string, std::vector<range>::iterator>::iterator lastfoundchrmit = lastfound.find(chrm);
	      if(lastfoundchrmit == lastfound.end()){ // first entry to this chrm
		      lastfound[chrm] = chrmit->second.begin();// start from begining
	      }
	      // set the hint to where the previous found is
	      searchhint = lastfound[chrm];
	      std::vector<range>::iterator itsearch = searchRange(searchhint, itend, rg);
	      // if found
	      if(itsearch != itend){// if found
		searchhint = itsearch;
		resultmap[tag].n_exreadsExcluded +=1;
		lastfound[chrm] = searchhint; // update search hint
		continue;
	      }
	    }

	  }
	}

	resultmap[tag].numTotal +=1;

	if(record1.IsMapped()) {
	    resultmap[tag].numMapped += 1;
	}

	if(record1.IsDuplicate()) {
	    resultmap[tag].numDuplicates +=1;
	}

	double gc = calcGC(record1.QueryBases);
	int ptn_count = countMotif(record1.QueryBases, opt::PATTERN, opt::PATTERN_REV);
	// when the read length exceeds 100bp, number of patterns might exceed the boundary
	if (ptn_count > ScanParameters::TEL_MOTIF_N-1){
	    continue;
	}
	resultmap[tag].telcounts[ptn_count]+=1;


	if(gc >= ScanParameters::GC_LOWERBOUND && gc <= ScanParameters::GC_UPPERBOUND){
	  // get index for GC bin.
	  int idx = floor((gc-ScanParameters::GC_LOWERBOUND)/ScanParameters::GC_BINSIZE);
	  assert(idx >=0 && idx <= ScanParameters::GC_BIN_N-1);
	  if(idx > ScanParameters::GC_BIN_N-1){
	    std::cerr << nprocessed << " GC:{"<< gc << "} telcounts:{"<< ptn_count <<"} GC bin index out of bound:" << idx << "\n";
	    exit(EXIT_FAILURE);
	  }
	  resultmap[tag].gccounts[idx]+=1;
	}

	nprocessed++;

	if( nprocessed%10000000 == 0){
	  std::cerr << "[scan] processed " << nprocessed << " reads \n" ;
	}
      }

      pBamReader->Close();
      delete pBamReader;

      // consider each BAM separately
      resultlist.push_back(resultmap);

      std::cerr << "[scan] total reads in BAM scanned " << ntotal << std::endl;
      std::cerr << "Completed scanning BAM\n";
  }

  if(opt::onebam){
    merge_results_by_readgroup(resultlist);
  }

  outputresults(resultlist);

  if(isExome){
    printlog(resultlist);
  }

  std::cerr << "Completed writing results\n";

  return 0;
}

void printlog(std::vector< std::map<std::string, ScanResults> > resultlist){

	for(size_t i =0; i< resultlist.size(); i++){
		auto rmap = resultlist[i];
		for(std::map<std::string, ScanResults>::iterator it= rmap.begin();
				it != rmap.end(); ++it){

			std::string rg = it ->first;
			ScanResults result = it -> second;
			std::cout << "BAM:" << rg << std::endl;
			std::cout << "	chr ID unmatched reads: " << result.n_exreadsChrUnmatched << std::endl;
			std::cout << "	exome reads excluded: " << result.n_exreadsExcluded << std::endl;
		}
	}

}


void printout(std::string rg, ScanResults result, std::ostream* pWriter){

	*pWriter << rg << ScanParameters::FIELD_SEP;
	*pWriter << result.lib << ScanParameters::FIELD_SEP;
	*pWriter << result.sample << ScanParameters::FIELD_SEP;
	*pWriter << result.numTotal << ScanParameters::FIELD_SEP;
	*pWriter << result.numMapped << ScanParameters::FIELD_SEP;
	*pWriter << result.numDuplicates << ScanParameters::FIELD_SEP;

	result.telLenEstimate = calcTelLength(result);
	if(result.telLenEstimate==-1){
		std::cerr << "Telomere length estimate unknown. No read was found with telomeric GC composition.\n";
		*pWriter << opt::unknown << ScanParameters::FIELD_SEP;
	}else if(result.telLenEstimate>1000000){
		std::cerr << "Telomere length estimate unknown. Possibly due to not enough representation of genome.\n";
		*pWriter << opt::unknown << ScanParameters::FIELD_SEP;
	}else if(result.telLenEstimate==0){
		std::cerr << "Telomere length estimate unknown. No read contains more than " << opt::tel_k << " telomere repeats.\n";
		*pWriter << opt::unknown << ScanParameters::FIELD_SEP;
	}
	else{
		*pWriter << result.telLenEstimate << ScanParameters::FIELD_SEP;
	}

	for (std::size_t j = 0, max = result.telcounts.size(); j != max; ++j){
		*pWriter << result.telcounts[j] << ScanParameters::FIELD_SEP;
	}
	for (std::size_t k = 0, max = result.gccounts.size(); k != max; ++k){
		*pWriter << result.gccounts[k] << ScanParameters::FIELD_SEP;
	}
	*pWriter << "\n";
}


int outputresults(std::vector< std::map<std::string, ScanResults> > resultlist){

	std::ostream* pWriter;
	bool tostdout = opt::outputfile.empty() ? true:false;
	if(tostdout){
		pWriter = &std::cout;
	}else{
		pWriter = createWriter(opt::outputfile);
	}

	if(opt::writerheader){
		Headers hd;
		for(size_t h=0; h<hd.headers.size();h++){
			*pWriter << hd.headers[h] << ScanParameters::FIELD_SEP;
		}
		*pWriter << "\n";
	}

	ScanResults mergedrs;
	std::string grpnames = "";

	for(size_t i=0; i < resultlist.size();++i){

		std::map<std::string, ScanResults> resultmap = resultlist[i];

		// if merge read groups, take weighted average for all measures
		bool domg = opt::mergerg && resultmap.size() > 1? true:false;

		for(std::map<std::string, ScanResults>::iterator it= resultmap.begin();
				it != resultmap.end(); ++it){

			std::string rg = it ->first;
			ScanResults result = it -> second;

			if(domg){
				if(grpnames.size()==0){
					grpnames += rg;
				}else{
					grpnames += "|"+rg;
				}

				mergedrs.sample = result.sample;
				mergedrs.numTotal += result.numTotal;
				mergedrs.numMapped += result.numMapped * result.numTotal;
				mergedrs.numDuplicates += result.numDuplicates * result.numTotal;
				mergedrs.telLenEstimate += calcTelLength(result)* result.numTotal;

				for (std::size_t j = 0, max = result.telcounts.size(); j != max; ++j){
					mergedrs.telcounts[j] += result.telcounts[j]* result.numTotal;
				}
				for (std::size_t k = 0, max = result.gccounts.size(); k != max; ++k){
					mergedrs.gccounts[k] += result.gccounts[k]* result.numTotal;
				}
				continue;
			}else{
				printout(rg, result, pWriter);
			}
		}

		//in this case calculate weighted average
		if(domg){

			mergedrs.numMapped /= mergedrs.numTotal;
			mergedrs.numDuplicates /= mergedrs.numTotal;
			mergedrs.telLenEstimate /= mergedrs.numTotal;

			for (std::size_t j = 0, max = mergedrs.telcounts.size(); j != max; ++j){
				mergedrs.telcounts[j] /= mergedrs.numTotal;
			}
			for (std::size_t k = 0, max = mergedrs.gccounts.size(); k != max; ++k){
				mergedrs.gccounts[k] /= mergedrs.numTotal;
			}

			mergedrs.numTotal  /= resultmap.size();

			printout(grpnames, mergedrs, pWriter);
		};
	}

	if(!tostdout){
		delete pWriter;
	}
	return 0;
}

double calcTelLength(ScanResults results){

	float acc = 0;
	for(std::size_t i=opt::tel_k, max = results.telcounts.size(); i !=max; ++i){
		acc += results.telcounts[i];
	}

	float gc_tel = 0;
	for(std::size_t i=0, max = results.gccounts.size(); i !=max; ++i){
		float gc1=ScanParameters::GC_LOWERBOUND + ScanParameters::GC_BINSIZE*i;
		float gc2=ScanParameters::GC_LOWERBOUND + ScanParameters::GC_BINSIZE*(i+1);
		if(gc1 >= ScanParameters::GC_TELOMERIC_LOWERBOUND && gc2 <= ScanParameters::GC_TELOMERIC_UPPERBOUND ){
			gc_tel += results.gccounts[i];
		}
	}

	if(gc_tel == 0){
		return -1;
	}

	return (acc/gc_tel)*float(ScanParameters::GENOME_LENGTH_AT_TEL_GC)/ScanParameters::LENGTH_UNIT/ScanParameters::TELOMERE_ENDS;
}



int countMotif(std::string &read, std::string pattern, std::string pattern_revcomp){

	int motifcount1 = 0;
	int motifcount2 = 0;

	size_t p1 = read.find(pattern, 0);
	while(p1 != std::string::npos)
	{
	    p1 = read.find(pattern,p1+pattern.size());
	    motifcount1 += 1;
	}

	size_t p2 = read.find(pattern_revcomp, 0);
	while(p2 != std::string::npos)
	{
	    p2 = read.find(pattern_revcomp,p2+pattern_revcomp.size());
	    motifcount2 += 1;
	}

	return motifcount1 > motifcount2? motifcount1:motifcount2;

}

double calcGC(const std::string& seq)
{
    double num_gc = 0.0f;
    double num_total = 0.0f;
    for(size_t i = 0; i < seq.size(); ++i)
    {
        if(seq[i] == 'C' || seq[i] == 'G')
            ++num_gc;
        ++num_total;
    }
    return num_gc / num_total;
}

void update_pattern(){

	opt::PATTERN = ScanParameters::PATTERN;
	opt::PATTERN_REV = ScanParameters::PATTERN_REVCOMP;

    // update total motif counts when read length and/or pattern have been specified by user
    ScanParameters::TEL_MOTIF_N = ScanParameters::READ_LENGTH/ScanParameters::PATTERN.size() +1;
	if(opt::tel_k < 1 or opt::tel_k > ScanParameters::TEL_MOTIF_N-1){
		std::cerr << "k out of bound. k must be an integer from 1 to " <<  ScanParameters::TEL_MOTIF_N-1 << "\n";
		exit(EXIT_FAILURE);
	}

}


//
// Handle command line arguments
//
void parseScanOptions(int argc, char** argv)
{

	std::string bamlistfile =  "";
	std::string rev = "";

	Headers hd;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'f':
            	arg >> bamlistfile; break;
            case 'o':
            	arg >> opt::outputfile; break;
            case 'H':
            	opt::writerheader=false; break;
            case 'm':
                opt::mergerg = true; break;
            case 'u':
                opt::ignorerg = true; break;
            case 'w':
                opt::onebam = true; break;
            case 'h':
        		for(size_t h=0; h<hd.headers.size();h++){
        			std::cout << hd.headers[h] << ScanParameters::FIELD_SEP;
        		}
        		std::cout << "\n";
        		exit(EXIT_SUCCESS);
            case 'k':
            	arg >> opt::tel_k;
            	break;
            case 'r':
				arg >> ScanParameters::READ_LENGTH;
				if(ScanParameters::READ_LENGTH <= 0 || ScanParameters::READ_LENGTH > 100000){
					std::cerr << "please specify valid read length that is greater than 0 and length than 100kb" << "\n";
					exit(EXIT_FAILURE);
				}
				break;
            case 't':
                {
                    int requestedThreads = 0;
                    arg >> requestedThreads;
                    if(!arg || requestedThreads < 1 || requestedThreads > 1024){
                        std::cerr << "threads must be an integer from 1 to 1024\n";
                        exit(EXIT_FAILURE);
                    }
                    opt::threads = static_cast<unsigned int>(requestedThreads);
                }
                break;
            case 'p':

				break;
            case 'z':
            	arg >> ScanParameters::PATTERN;
				ScanParameters::PATTERN_REVCOMP = reverseComplement(ScanParameters::PATTERN);
				std::cerr << "use user specified pattern " <<  ScanParameters::PATTERN << "\n";
				std::cerr << "reverse complement " <<  ScanParameters::PATTERN_REVCOMP << "\n";
            	break;
            case 'e':
            	arg >> opt::exomebedfile;
            	opt::exomebed = readBedAsVector(opt::exomebedfile);
            	std::cout << "loaded "<< opt::exomebed.size() << " exome regions \n"<< std::endl;
//            	std::cout << opt::exomebed << "\n";
            	break;

            case OPT_HELP:
                std::cout << TELSEQ_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << TELSEQ_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    update_pattern();

    // deal with cases of API usage:
    // telseq a.bam b.bam c.bam ...
    // | telseq
    // telseq -f

    if (argc - optind < 1) // no argument specified
    {
    	// check if it is from pipe
		if(!isatty(fileno(stdin))){
			std::string line;
			while (std::getline(std::cin, line))
			{
				if(line.empty()){
					continue;
				}
				opt::bamlist.push_back(line);
			}
		}else if(bamlistfile.empty() ){ // check if not from a pipe, -f must be spceified
//    		std::cerr << SUBPROGRAM ": No BAM specified. Please specify BAM either directly, by using -f option or piping BAM file path.\n";
    		std::cout << TELSEQ_USAGE_MESSAGE;
    		exit(EXIT_FAILURE);
    	}
    }

    else if (argc - optind >= 1) // if arguments are specified
    {
    	// -f has higher priority, when specified, ignore arguments.
    	if(bamlistfile.empty()){
    		for(int i = optind; i < argc; ++i ){
    		    opt::bamlist.push_back(argv[i]);
    		}
    	}
    }

    // read in bamlist
    if(!bamlistfile.empty()){
        size_t filesize = getFilesize(bamlistfile);
        if(filesize == 0)
        {
            std::cerr << PROGRAM_BIN ": BAMLIST file specified by -f is empty\n";
            exit(EXIT_FAILURE);
        }

        std::istream* pReader = createReader(bamlistfile);
        std::string line;

        while(getline(*pReader, line))
        {
        	if(line.empty()){
        		continue;
        	}
            opt::bamlist.push_back(line);
        }

        size_t bamsize = opt::bamlist.size();
        if(bamsize == 0 ){
            std::cerr << PROGRAM_BIN ": Could not find any sample in BAMLIST file specified.\n";
            exit(EXIT_FAILURE);
        }
        delete pReader;
    }
}



int main(int argc, char** argv)
{
	Timer* pTimer = new Timer("scan BAM");
	parseScanOptions(argc, argv);
    const int status = scanBam();
    delete pTimer;
    return status;
}
