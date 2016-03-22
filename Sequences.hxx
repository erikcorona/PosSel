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

namespace Gen
{
    class Sequence : public std::string
    {
    public:
        Sequence(std::string s) : std::string(s){};

    private:

    };

    class Sequences
    {
        friend class HaploidSequences;
        friend class DiploidSequences;
    public:

        /**
         * Adds one sequence to the list of sequences
         * @param seq the haplotype sequence
         * @param orgID the ID of the organism
         */
        void addSeq(std::string seq)
        {
            seqs.push_back(std::make_shared<Sequence>(seq));
        }

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
                this->addSeq(*seq.get());
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

    class DiploidSequences : public Sequences {
    public:

        DiploidSequences(){}
        DiploidSequences(DiploidSequences& other)
        {
            for(auto& seq : other.seqs)
                addSeq(*seq.get());
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
