//
// Created by ecorona on 3/21/16.
//

#ifndef POSSEL_SEQUENCES_HXX
#define POSSEL_SEQUENCES_HXX

#include <cmath>
#include <boost/random.hpp>
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

        Sequence(){}
        Sequence(std::string s)
        {
            seq = s;
        };

        void addAllele(char allele)
        {
            if(allele != 'A' && allele != 'C' && allele != 'G' && allele != 'T')
                assert(false);
            seq += allele;
        }

        virtual std::string getString()=0;
        virtual char& operator[](std::size_t i)=0;
        virtual std::size_t size()=0;
    protected:
        std::string seq;
    };

    class StrSequence : public Sequence{
    public:
        StrSequence(std::string s) : Sequence(s) { }
        char& operator[](std::size_t i) { return seq[i];}
        std::size_t size(){return seq.size();}
        std::string getString(){return seq;}
    };

    class LongSequence : public Sequence
    {
    public:
        LongSequence(std::vector<std::string>& rsids, std::vector<int>& positions)
        : ids{rsids}, pos{positions}
        {
            startI = 0;
            endI   = 0; // exclusive
        }

        std::string getString(){
            return seq.substr(static_cast<std::size_t>(startI), static_cast<std::size_t>(endI-startI));
        }

        std::size_t size(){return static_cast<std::size_t>(endI - startI);}
        auto trueSize(){return seq.size();}
        char& operator[](std::size_t i)
        {
            return seq[startI+i];
        }

        /**
         * @param startIndex inclusive
         * @param endIndex exclusive
         */
        void setWindowByIndex(int startIndex, int endIndex)
        {
            startI = startIndex;
            endI = endIndex;
        }
    private:
        int startI, endI;
        std::vector<std::string>& ids;
        std::vector<int>& pos;
    };

    class Sequences
    {
        friend class HaploidSequences;
        friend class DiploidSequences;
    public:

        char allele(std::size_t seqIndex, std::size_t alleleIndex)
        {
            return (*seqs[seqIndex])[alleleIndex];
        }

        void showSequences()
        {
            for(auto& s : seqs)
            {
                for(std::size_t i{0}; i < seqLength(); i++)
                    std::cout << (*s)[i];
                std::cout << std::endl;
            }
        }

        std::string getString(std::size_t i)
        {
            return seqs[i]->getString();
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
                if(seqs[0]->size() != sequence->size())
                    return false;
            return true;
        }

        std::size_t seqLength(){return seqs[0]->size();}

        template<typename Fxn>
        std::vector<std::unique_ptr<Sequences>> sortedSamples(std::size_t n, Fxn sorter)
        {
            std::vector<std::unique_ptr<Sequences>> vec(n);

            for(auto& v : vec)
                v = sample();
            std::sort(vec.begin(), vec.end(), sorter);
            return vec;
        }

        bool isSegregating(std::size_t i)
        {
            char ref = (*seqs[0])[i];
            for(auto& seq : seqs)
            {
                if ((*seq)[i] != ref)
                    return true;
            }
            return false;
        }

        auto nSequences(){ return seqs.size();}
        virtual std::size_t nOrganisms()=0;
        virtual void addOrganism(Sequences& sequences, int orgI)=0;

    protected:
        std::vector<std::shared_ptr<Sequence>> seqs;

        virtual std::unique_ptr<Sequences> uniq_new()=0;

        std::unique_ptr<Sequences> randomSitesClone()
        {
            static boost::random::mt19937 rng(2);
            boost::random::uniform_int_distribution<> sites(0,static_cast<int>(seqLength())-1);
            auto copy = clone();
            for(std::size_t i = 0; i < seqLength(); i++)
            {
                int siteIndex = sites(rng);
                for(std::size_t seqI = 0; seqI < seqs.size(); seqI++)
                    (*copy->seqs[seqI].get())[i] = (*seqs[seqI].get())[siteIndex];
            }
            return copy;
        }
        std::unique_ptr<Sequences> sample()
        {
            static boost::random::mt19937 rng(1);
            boost::random::uniform_int_distribution<> org(0,static_cast<int>(nOrganisms())-1);

            auto sequences = uniq_new();

            auto tmp = randomSitesClone();
            for(std::size_t i = 0; i < nOrganisms(); i++)
                sequences->addOrganism(*tmp, org(rng));
            return sequences;
        }

        virtual std::unique_ptr<Sequences> clone()=0;
    };

    class HaploidSequences : public Sequences
    {
    public:

        HaploidSequences(){}
        HaploidSequences(HaploidSequences& other)
        {
            for(auto& seq : other.seqs)
                addSeq(seq);
        }

        void addOrganism(Sequences& sequences, int i)
        {
            seqs.push_back(sequences.seqs[i]);
        }

        std::size_t nOrganisms() { return seqs.size();}

    protected:
        std::unique_ptr<Sequences> clone()
        {
            return std::unique_ptr<Sequences>(new HaploidSequences(*this));
        }

    private:
        std::unique_ptr<Sequences> uniq_new()
        {
            return std::unique_ptr<Sequences>(new HaploidSequences());
        }
    };

    class PositionSequences
    {

    public:
        void setBuild(std::string theBuild)
        {
            build = theBuild;
        }

        void setChr(short chromosome){ chr = chromosome;}

    protected:
        short chr;
        std::vector<int> pos;
        std::string build;
    };

    class GeneticMap : public PositionSequences
    {
    public:
        GeneticMap()
        {
            ensureDataExists();
//            readInData();
        }

        double getcM(int chrPos)
        {
            // binary search for ---|*|----
            int i = win::getIndexBehindOrEqual(pos,chrPos);
            if(chrPos == pos[i])
                return cM[i];
            assert(i+1 < pos.size());
            return cM[i] + (chrPos-pos[i])*(cM[i+1] - cM[i])/(pos[i+1] - pos[i]); // extrapolate cM
        }

        double distance(int start, int end)
        {
            assert(end > start);
            double val = getcM(start) - getcM(end);
            return val >= 0 ? val : -val;
        }

    private:

        std::vector<double> cM; // first index is chr, 2nd is the position
        void ensureDataExists()
        {
            if(!dataExists(dataSetPath()))
                download(dataSetPath());
        }

        void download(std::string path)
        {
            system(std::string("wget ") +
                           std::string("https://hapmap.ncbi.nlm.nih.gov/downloads/recombination/latest/rates/genetic_map_chr") +
                           std::to_string(chr) +
                           std::string("_b36.txt") +
                           std::string("-O ") +
                           path
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



class HapMapSequences : public HaploidSequences, public PositionSequences
    {
    public:

        HapMapSequences(std::string fileName)
        {
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
                    seqs[i] = std::shared_ptr<Sequence>(new LongSequence(rsids, pos));

                while(!file.eof())
                {
                    getline(file,line);

                    if(line.size() < 1) // skip blank lines
                        continue;
                    boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
                    std::vector<std::string> row(tokens.begin(), tokens.end());

                    assert(row.size() > 0);
                    rsids.push_back(row[0]);
                    pos.push_back(atoi(row[1].c_str()));

                    for(std::size_t i = 2; i < row.size(); i++)
                        seqs[i-2]->addAllele(row[i][0]);
                }
            }
        }

        auto trueSeqLength()
        {
            return static_cast<LongSequence*>(seqs[0].get())->trueSize();
        }

        /**
         * @param startIndex inclusive
         * @param endIndex exclusive
         */
        template<typename WholeNumber>
        void setWindowByIndex(WholeNumber startIndex, WholeNumber endIndex)
        {
            for(auto& seq : seqs)
            {
                LongSequence* ls;
                ls = static_cast<LongSequence*>(seq.get());
                ls->setWindowByIndex(startIndex, endIndex);
            }
        }

        void setWindowByPos(int startPos, int endPos)
        {
            int startI = win::getIndexBehindOrEqual(pos,startPos);
            int endI = win::getIndexInfront(pos,endPos);
            std::cerr << "s:" << startI << "\t e:" << endI << std::endl;
            assert(startI >= 0 && endI >= 0 && endI != startI);
            setWindowByIndex(startI, endI);
        }

    private:
        std::vector<std::string> rsids;
    };

    class DiploidSequences : public Sequences
    {
    public:

        DiploidSequences(){}
        DiploidSequences(DiploidSequences& other)
        {
            for(auto& seq : other.seqs)
                addSeq(seq);
        }

        void addOrganism(Sequences &sequences, int i)
        {
            seqs.push_back(sequences.seqs[i*2]);
            seqs.push_back(sequences.seqs[i*2 + 1]);
        }

    protected:
        std::unique_ptr<Sequences> clone()
        {
            return std::unique_ptr<Sequences>(new DiploidSequences(*this));
        }

    private:
        std::unique_ptr<Sequences> uniq_new()
        {
            return std::unique_ptr<Sequences>(new DiploidSequences());
        }

        std::size_t nOrganisms() { return seqs.size()/2;}
    };
}

#endif //POSSEL_SEQUENCES_HXX
