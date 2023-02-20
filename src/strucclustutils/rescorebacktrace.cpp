#include "DBReader.h"
#include "IndexReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Alignment.h"
#include "structureto3diseqdist.h"
#include "StructureSmithWaterman.h"
#include "StructureUtil.h"
#include "TMaligner.h"
#include "Coordinate16.h"
#include "LDDT.h"

#ifdef OPENMP
#include <omp.h>
#endif


int rescorebacktrace(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader qdbrAA(par.db1, par.threads, IndexReader::SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    IndexReader qdbr3Di(StructureUtil::getIndexWithSuffix(par.db1, "_ss"), par.threads, IndexReader::SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);

    IndexReader *t3DiDbr = NULL;
    IndexReader *tAADbr = NULL;
    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        t3DiDbr = &qdbr3Di;
        tAADbr = &qdbrAA;
    } else {
        tAADbr = new IndexReader(par.db2, par.threads, IndexReader::SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
        t3DiDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db2, "_ss"), par.threads, IndexReader::SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    }

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed,  Parameters::DBTYPE_ALIGNMENT_RES);
    dbw.open();


    SubstitutionMatrix subMat3Di(par.scoringMatrixFile.values.aminoacid().c_str(), 2.1, par.scoreBias);
    std::string blosum;
    for (size_t i = 0; i < par.substitutionMatrices.size(); i++) {
        if (par.substitutionMatrices[i].name == "blosum62.out") {
            std::string matrixData((const char *)par.substitutionMatrices[i].subMatData, par.substitutionMatrices[i].subMatDataLen);
            std::string matrixName = par.substitutionMatrices[i].name;
            char * serializedMatrix = BaseMatrix::serialize(matrixName, matrixData);
            blosum.assign(serializedMatrix);
            free(serializedMatrix);
            break;
        }
    }
    SubstitutionMatrix subMatAA(blosum.c_str(), 1.4, par.scoreBias);
    //temporary output file
    Debug::Progress progress(resultReader.getSize());


#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<Matcher::result_t> alignmentResult;
        TMaligner *tmaligner = NULL;
        Sequence qSeqAA(par.maxSeqLen, qdbrAA.getDbtype(), (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence qSeq3Di(par.maxSeqLen, qdbr3Di.getDbtype(), (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        Sequence tSeqAA(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence tSeq3Di(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        std::string backtrace;
        char buffer[1024+32768];
        std::string resultBuffer;

        std::vector<Matcher::result_t> alignments;

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();
            char *data = resultReader.getData(id, thread_idx);

            Matcher::readAlignmentResults(alignments, data, false);
            size_t queryKey = resultReader.getDbKey(id);

            unsigned int queryId3Di = qdbr3Di.sequenceReader->getId(queryKey);
            unsigned int queryIdAA = qdbrAA.sequenceReader->getId(queryKey);
            char *querySeqAA = qdbrAA.sequenceReader->getData(queryIdAA, thread_idx);
            char *querySeq3Di = qdbr3Di.sequenceReader->getData(queryId3Di, thread_idx);
            unsigned int querySeqLen3Di = qdbr3Di.sequenceReader->getSeqLen(queryId3Di);
            unsigned int querySeqLenAA = qdbrAA.sequenceReader->getSeqLen(queryIdAA);

            qSeq3Di.mapSequence(queryId3Di, queryKey, querySeq3Di, querySeqLen3Di);
            qSeqAA.mapSequence(queryIdAA, queryKey, querySeqAA, querySeqLenAA);

            int gapOpen, gapExtend;
            unsigned int gapIcount, gapDcount;
            int rescore, gapScore;
            gapOpen = par.gapOpen.values.aminoacid();
            gapExtend = par.gapExtend.values.aminoacid();

            for(size_t alnIdx = 0; alnIdx < alignments.size(); alnIdx++){
                //data = Util::skipLine(data);
                const unsigned int dbKey = alignments[alnIdx].dbKey;
                unsigned int targetId3Di = t3DiDbr->sequenceReader->getId(dbKey);
                unsigned int targetIdAA = tAADbr->sequenceReader->getId(dbKey);
                //const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB))? true : false;

                char * targetSeq3Di = t3DiDbr->sequenceReader->getData(targetId3Di, thread_idx);
                char * targetSeqAA = tAADbr->sequenceReader->getData(targetIdAA, thread_idx);
                const int targetSeqLen3Di = static_cast<int>(t3DiDbr->sequenceReader->getSeqLen(targetId3Di));
                const int targetSeqLenAA = static_cast<int>(tAADbr->sequenceReader->getSeqLen(targetIdAA));

                tSeq3Di.mapSequence(targetId3Di, dbKey, targetSeq3Di, targetSeqLen3Di);
                tSeqAA.mapSequence(targetIdAA, dbKey, targetSeqAA, targetSeqLenAA);

                rescore = 0;

                int qpos = isdigit(alignments[alnIdx].qStartPos);
                int dbpos = isdigit(alignments[alnIdx].dbStartPos);

                for(size_t  j = 0; j < alignments[alnIdx].backtrace.size(); j++){

                    //if match add 3Di and AA score as rescore value
                    char qAALetter = qSeqAA.numSequence[qpos];
                    char dbAALetter =  tSeqAA.numSequence[dbpos];
                    char q3DiLetter = qSeq3Di.numSequence[qpos];
                    char db3DiLetter =  tSeq3Di.numSequence[dbpos];

                    switch (alignments[alnIdx].backtrace[j]) {
                        case 'M':
                            rescore += subMatAA.subMatrix[qAALetter][dbAALetter] + subMat3Di.subMatrix[q3DiLetter][db3DiLetter];
                            gapIcount = 0; gapDcount = 0;
                            qpos++;
                            dbpos++;
                            break;
                        case 'D':
                            gapIcount = 0;
                            gapScore = (gapDcount == 0) ? gapOpen : gapExtend;
                            rescore += gapScore;
                            gapDcount += 1;
                            qpos++;
                            break;
                        case 'I':
                            gapDcount = 0;
                            gapScore = (gapIcount == 0) ? gapOpen : gapExtend;
                            rescore += gapScore;
                            gapIcount += 1;
                            dbpos++;
                            break;
                    }
                }

                alignments[alnIdx].score = rescore;
                size_t len = Matcher::resultToBuffer(buffer, alignments[alnIdx], par.addBacktrace);
                resultBuffer.append(buffer, len);

                dbw.writeData(resultBuffer.c_str(), resultBuffer.length(), queryKey, thread_idx);
                resultBuffer.clear();
                alignmentResult.clear();
            }
        }
    }


    dbw.close();
    resultReader.close();


    if (sameDB == false) {
        delete t3DiDbr;
        delete tAADbr;
    }

    return EXIT_SUCCESS;
}
