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

        std::size_t size(){return static_cast<std::size_t>(endI - startI);}
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

        void showSequences()
        {
            for(auto& s : seqs)
            {
                for(std::size_t i{0}; i < seqLength(); i++)
                    std::cout << (*s)[i];
                std::cout << std::endl;
            }
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

        double tajPI()
        {
            int diff{0};
            for(std::size_t i = 0; i < seqs.size(); i++)
                for(std::size_t k = i+1; k < seqs.size(); k++)
                    for(std::size_t allele{0}; allele < seqLength(); allele++)
                        diff += (*seqs[i])[allele] != (*seqs[k])[allele];

            return static_cast<double>(diff)/(seqs.size()*(seqs.size()-1.0)/2);
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

        int tajS()
        {
            int segs{0};
            for(std::size_t i = 0; i < seqLength(); i++)
                segs += isSegregating(i);
            return segs;
        }

        double tajD()
        {
            double pi = tajPI();
            double s = tajS();
            return (pi-s/a1())/(std::sqrt(e1()*s+e2()*s*(s-1)));
        }

        virtual std::size_t nOrganisms()=0;
        virtual void addOrganism(Sequences& sequences, int orgI)=0;

        void print()
        {
            std::cout << seqs.size() << " sequences, " << seqLength() << " sites" << std::endl;
            std::cout << "pi: " << tajPI() << std::endl;
            std::cout << "Segregating sites: " << tajS() << "/" << seqLength() << std::endl;
            std::cout << "theta_hat[estimated from S]: " << tajS()/a1() << std::endl;
            std::cout << "Tajima's D: " << tajD() << std::endl;
            std::cout << "a1=" << a1() << " a2=" << a2() << " b1=" << b1() << " b2=" << b2() << std::endl;
            std::cout << "c1=" << c1() << " c2=" << c2() << " e1=" << e1() << " e2=" << e2() << std::endl;
        }

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
    private:

        // a1 is described in tajmad1.pdf
        double a1()
        {
            double sum{0};
            for(std::size_t i = 1; i < seqs.size(); i++)
                sum += 1.0/i;
            return sum;
        }

        double a2()
        {
            double sum{0};
            for(std::size_t i = 1; i < seqs.size(); i++)
                sum += 1.0/(i*i);
            return sum;
        }

        double b1()
        {
            double n = static_cast<double>(seqs.size());
            return (n+1.0)/(3.0*(n-1));
        }

        double b2()
        {
            double n = static_cast<double>(seqs.size());
            return (2*(n*n+n+3))/(9*n*(n-1));
        }

        double c1() {            return b1()-1.0/a1();        }

        double c2()
        {
            double n = static_cast<double>(seqs.size());
            return b2() - (n+2)/(a1()*n) + a2()/(a1()*a1());
        }

        double e1()        {            return c1()/a1();        }
        double e2()        {            return c2()/(a1()*a1()+a2());        }
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

    class HapMapSequences : public HaploidSequences
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
            int startI{-2}, endI{-2};
            for(std::size_t i = 0; i < pos.size(); i++)
            {
                if (pos[i] == startPos) {
                    startI = static_cast<int>(i);
                    break;
                }

                if(pos[i] > startPos)
                {
                    startI = static_cast<int>(i-1);
                    break;
                }
            }

            if(startI < 0)
                startI = 0;

            for(std::size_t i = 0; i < pos.size(); i++)
            {
                if(pos[i] == endPos)
                {
                    endI = static_cast<int>(i + 1);
                    break;
                }

                if(pos[i] > endPos)
                {
                    endI = static_cast<int>(i);
                    break;
                }
            }

            if(endPos > pos.back())
                endI = static_cast<int>(pos.size());
            std::cerr << "s:" << startI << "\t e:" << endI << std::endl;
            assert(startI >= 0 && endI >= 0 && endI != startI);
            setWindowByIndex(startI, endI);
            assert(seqLength() > 0);
        }
    private:
        std::vector<std::string> rsids;
        std::vector<int> pos;
    };

    class DiploidSequences : public Sequences {
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
