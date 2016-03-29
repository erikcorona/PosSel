//
// Created by ecorona on 3/21/16.
//

#ifndef POSSEL_SEQUENCES_HXX
#define POSSEL_SEQUENCES_HXX

#include <cmath>
#include <boost/random.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <iostream>
#include <string>
#include <unordered_map>
#include <memory>
#include <fstream>
#include <boost/tokenizer.hpp>
#include "WindowAlgos.hxx"

//@todo bootstrap the value (PI - S/a1) to get an empirical pr(diff) btw PI & S/a1
namespace Gen
{
    class Sequence
    {
    public:
        Sequence(std::string s)
        {
            seq = s;
//            startI = 0;
//            endI = static_cast<int>(s.size());
        }

        Sequence()
        {
//            startI = 0;
//            endI   = 0; // exclusive
        }

        bool operator==(Sequence* b)
        {
//            assert(startI == b->startI && endI == b->endI);
//            if(startI != b->startI)
//                return false;
//
//            if(endI != b->endI)
//                return false;

            return seq == b->seq;
        }

        bool operator !=(Sequence* b)
        {
            return !(this->operator==(b));
        }

        void addAllele(char allele)
        {
            if(allele != 'A' && allele != 'C' && allele != 'G' && allele != 'T')
                assert(false);
            seq += allele;
        }

        std::string getString(std::size_t startI, std::size_t endI){
            return seq.substr(static_cast<std::size_t>(startI), static_cast<std::size_t>(endI-startI));
        }

//        std::size_t size(){return static_cast<std::size_t>(endI - startI);}
        auto trueSize(){return seq.size();}

        template<typename Index>
        char& allele(Index i, std::size_t startOffset)
        {
            return seq[startOffset+i];
        }

        /**
         * @param startIndex inclusive
         * @param endIndex exclusive
         */
//        void setWindowByIndex(int startIndex, int endIndex)
//        {
//            startI = startIndex;
//            endI = endIndex;
//        }

        /**
         * @param start inclusive
         * @param end exclusive
         */
//        void subSeqs(std::size_t start, std::size_t end)
//        {
//            setWindowByIndex(start,end);
//        }

    private:
//        int startI, endI;
        std::string seq;
    };

    class PositionSequences
    {

    public:
        PositionSequences(PositionSequences* other)
        {
            pos = other->pos;
            chr = other->chr;
            build = other->build;
        }

        bool operator==(PositionSequences* b)
        {
            assert(pos->size() == b->pos->size());
            if (pos->size() != b->pos->size())
                return false;

            assert(chr == b->chr);
            if(chr != b->chr)
                return false;

            assert(build == b->build);
            if(build != b->build)
                return false;

            for(std::size_t i = 0; i < pos->size(); i++)
            {

                assert((*pos)[i] == (*(b->pos))[i]);
                if ((*pos)[i] != (*(b->pos))[i])
                    return false;
            }
            return true;
        }

        PositionSequences()
        {
            pos = std::make_shared<std::vector<int>>();
        }

        void setBuild(std::string theBuild)
        {
            build = theBuild;
        }

        void setChr(short chromosome){ chr = chromosome;}
        auto getChr(){ return chr;}

    protected:
        short chr;
        std::shared_ptr<std::vector<int>> pos;
        std::string build;
    };

    class Sequences : public PositionSequences
    {
        friend class HaploidSequences;
        friend class DiploidSequences;
    public:

        Sequences()
        {
            seqs = std::vector<std::shared_ptr<Sequence>>();
            startI = 0;
            endI   = 0;
        }

        Sequences(Sequences* other)
        : PositionSequences(static_cast<PositionSequences*>(other))
        {
            startI = other->startI;
            endI   = other->endI;
            for(auto s : seqs)
                seqs.push_back(s);
        }

        bool operator==(Sequences* b)
        {
            assert(b->seqs.size() == seqs.size());
            if(b->seqs.size() != seqs.size())
                return false;

            if(startI != b->endI)
                return false;

            if(endI != b->endI)
                return false;

            for(std::size_t i = 0; i < seqs.size(); i++)
            {
                assert(seqs[i]->operator==(b->seqs[i].get()));
                if (seqs[i]->operator!=(b->seqs[i].get()))
                    return false;
            }

            assert(static_cast<PositionSequences*>(b)->operator==(static_cast<PositionSequences*>(this)));
            return static_cast<PositionSequences*>(b)->operator==(static_cast<PositionSequences*>(this));
        }

        char allele(std::size_t seqIndex, std::size_t alleleIndex)
        {
            return seqs[seqIndex]->allele(alleleIndex,startI);
        }

        void showSequences()
        {
            for(auto& s : seqs)
            {
                for(std::size_t i{0}; i < seqLength(); i++)
                    std::cout << s->allele(i, startI);
                std::cout << std::endl;
            }
        }

        auto getAllSeqs(){ return seqs;}
        void setSequences(std::vector<std::shared_ptr<Sequence>>& theSeqs)
        {
            seqs = theSeqs;
        }

        std::string getString(std::size_t seqIndex)
        {
            return seqs[seqIndex]->getString(startI, endI);
        }

        /**
         * Adds one sequence to the list of sequences
         * @param seq the haplotype sequence
         * @param orgID the ID of the organism
         */
        void addSeq(std::shared_ptr<Gen::Sequence> seq)
        {
            seqs.push_back(seq);
        }

        /**
         * Makes sure all sequences have the same length
         */
        bool isValid()
        {
            for(auto& sequence : seqs)
                if(sequence->trueSize() <= endI)
                    return false;
            return true;
        }

        std::size_t trueSeqLength(){return seqs[0]->trueSize();}
        std::size_t seqLength(){return endI - startI;}

        bool isSegregating(std::size_t i)
        {
            char ref = seqs[0]->allele(i, startI);
            for(auto& seq : seqs)
            {
                if (seq->allele(i,startI) != ref)
                    return true;
            }
            return false;
        }

        auto nSequences(){ return seqs.size();}

        template<typename Index>
        auto getPos(Index i) {
            return (*pos)[i];
        }

        template<typename Index>
        Sequences* setWindowByIndex(Index start, Index end)
        {
            startI = start;
            endI   = end;
            return this;
        }

        template<typename Fxn>
        void filter(Fxn f)
        {
            boost::remove_erase_if(seqs,f);
        }

        auto getStartI(){ return startI;}
        auto getEndI()  {return endI;}
        virtual std::size_t nOrganisms()=0;

        /**
         * Creates a shallow copy of the object. Primitives are fully copied but the sequences
         * are shared. It's best to treat sequences as read only to avoid problems  .
         */
        virtual std::shared_ptr<Sequences> clone()=0;

    protected:
        std::vector<std::shared_ptr<Sequence>> seqs;

        virtual std::shared_ptr<Sequences> uniq_new()=0;

        std::shared_ptr<Sequences> randomSitesClone()
        {
            static boost::random::mt19937 rng(2);
            boost::random::uniform_int_distribution<> sites(0,static_cast<int>(seqLength())-1);
            auto copy = clone();
            for(std::size_t i = 0; i < seqLength(); i++)
            {
                int siteIndex = sites(rng);
                for(std::size_t seqI = 0; seqI < seqs.size(); seqI++)
                    copy->seqs[seqI]->allele(i, startI) = seqs[seqI]->allele(siteIndex, startI);
            }
            return copy;
        }

    private:
        std::size_t startI, endI;

    };

    class HaploidSequences : public Sequences
    {
    public:

        HaploidSequences(){}
        HaploidSequences(HaploidSequences* other)
        : Sequences(static_cast<Sequences*>(other))
        { }

        void addOrganism(Sequences& sequences, int i)
        {
            seqs.push_back(sequences.seqs[i]);
        }

        std::size_t nOrganisms() { return seqs.size();}

        std::shared_ptr<Sequences> clone()
        {
            return std::shared_ptr<Sequences>(new HaploidSequences(*this));
        }

    private:
        std::shared_ptr<Sequences> uniq_new()
        {
            return std::shared_ptr<Sequences>(new HaploidSequences());
        }
    };

    class GeneticMap : public PositionSequences
    {
    public:
        GeneticMap(short chromosome)
        {
            chr = chromosome;
            ensureDataExists();
            readInData();
        }

        double getcM(int chrPos)
        {
            // binary search for ---|*|----
            int i = win::getIndexBehindOrEqual(*pos,chrPos);
            if(chrPos == (*pos)[i])
                return cM[i];
            assert(i+1 < static_cast<int>(pos->size()));
            return cM[i] + (chrPos-(*pos)[i])*(cM[i+1] - cM[i])/((*pos)[i+1] - (*pos)[i]); // extrapolate cM
        }

        double distance(int start, int end)
        {
            assert(end >= start);
            double val = getcM(start) - getcM(end);
            return val >= 0 ? val : -val;
        }

        bool dataExists(){ return dataExists(dataSetPath());}

    private:

        void readInData()
        {
            std::ifstream file;
            file.open(dataSetPath());
            if(file.is_open())
            {
                boost::char_separator<char> sep(" ");
                std::string line;
                if(!file.eof()) getline(file, line); // skip header

                while(!file.eof())
                {
                    getline(file,line);

                    if(line.size() < 1) // skip blank lines
                        continue;
                    boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
                    std::vector<std::string> row(tokens.begin(), tokens.end());
                    assert(row.size() > 0);
                    pos->push_back(atoi(row[0].c_str()));
                    cM.push_back(atof(row[2].c_str()));
                }
            }
        }

        std::vector<double> cM; // first index is chr, 2nd is the position
        void ensureDataExists()
        {
            if(!dataExists(dataSetPath()))
                download(dataSetPath());
        }

        int download(std::string path)
        {
            return system(std::string(std::string("wget ") +
                           std::string("https://hapmap.ncbi.nlm.nih.gov/downloads/recombination/latest/rates/genetic_map_chr") +
                           std::to_string(chr) +
                           std::string("_b36.txt ") +
                           std::string("-O ") +
                           path).c_str()
            );
        }

        bool dataExists(std::string path)
        {
            std::ifstream f(path.c_str());
            bool good = f.good();
            f.close();
            return good;
        }

        std::string dataSetPath()
        {
            using namespace std;
            return string("data/hapmapGeneticMap/genetic_map_chr") + to_string(chr) + string("_b36.txt");
        }
    };

    class DiploidSequences : public Sequences
    {
    public:

        DiploidSequences(){}
        DiploidSequences(DiploidSequences* other)
        : Sequences(static_cast<Sequences*>(other)){

        }

        bool operator==(DiploidSequences* b)
        {
            return static_cast<Sequences*>(b)->operator==(static_cast<Sequences*>(this));
        }

        std::shared_ptr<Sequences> clone()
        {
            return std::shared_ptr<Sequences>(new DiploidSequences(*this));
        }

    private:
        std::shared_ptr<Sequences> uniq_new()
        {
            return std::shared_ptr<Sequences>(new DiploidSequences());
        }

        std::size_t nOrganisms() { return seqs.size()/2;}
    };

class HapMapSequences : public DiploidSequences
    {
    public:

        HapMapSequences(HapMapSequences* toCopy)
        : DiploidSequences(static_cast<DiploidSequences*>(toCopy)){
            rsids = toCopy->rsids;
        }

        bool operator==(HapMapSequences* b)
        {
            for(std::size_t i = 0; i < b->rsids->size(); i++)
            {
                assert((*rsids)[i] == (*(b->rsids))[i]);
                if ((*rsids)[i] != (*(b->rsids))[i])
                    return false;
            }

            assert(static_cast<DiploidSequences*>(b)->operator==(static_cast<DiploidSequences*>(this)));
            return static_cast<DiploidSequences*>(b)->operator==(static_cast<DiploidSequences*>(this));
        }

        bool operator!=(HapMapSequences* b) {
            return !((*this) == b);
        }

//    HapMapSequences(HapMapSequences* toCopy)
//            : DiploidSequences(static_cast<DiploidSequences*>(toCopy)){
//        rsids = toCopy->rsids;
//        pos   = toCopy->pos;
//        std::cerr << this->getPos(0) << std::endl;
//        chr   = toCopy->chr;
//        build = toCopy->build;
//        for(auto& s : toCopy->seqs)
//            seqs.push_back(s);
//    }

        HapMapSequences(std::string fileName)
        {
            rsids = std::make_shared<std::vector<std::string>>();
            std::ifstream file;
            file.open(fileName);
            if(file.is_open())
            {
                boost::char_separator<char> sep(" ");
                std::string line;
                unsigned long numSeqs{0};
                if(!file.eof())
                {
                    getline(file, line); // get num individuals from header
                    boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
                    std::vector<std::string> toks(tokens.begin(), tokens.end());
                    numSeqs = toks.size()-2; // RSID and POS are not sequences
                }

                seqs.resize(numSeqs);
                for(std::size_t i = 0; i < numSeqs; i++)
                    seqs[i] = std::shared_ptr<Sequence>(new Sequence());

                while(!file.eof())
                {
                    getline(file,line);

                    if(line.size() < 1) // skip blank lines
                        continue;
                    boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
                    std::vector<std::string> row(tokens.begin(), tokens.end());

                    assert(row.size() > 0);
                    rsids->push_back(row[0]);
                    pos->push_back(atoi(row[1].c_str()));

                    for(std::size_t i = 2; i < row.size(); i++)
                        seqs[i-2]->addAllele(row[i][0]);
                }
            }
        }

        void setWindowByPos(int startPos, int endPos)
        {
            int startI = win::getIndexBehindOrEqual(*pos,startPos);
            int endI = win::getIndexInfront(*pos,endPos);
            std::cerr << "s:" << startI << "\t e:" << endI << std::endl;
            assert(startI >= 0 && endI >= 0 && endI != startI);
            setWindowByIndex(startI, endI);
        }

        std::shared_ptr<Sequences> clone()
        {
            return std::shared_ptr<Sequences>(new HapMapSequences(*this));
        }

    private:
        std::shared_ptr<std::vector<std::string>> rsids;
    };


}

#endif //POSSEL_SEQUENCES_HXX
