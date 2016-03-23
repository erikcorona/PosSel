//
// Created by ecorona on 3/23/16.
//

#ifndef POSSEL_SEQUENCESALGOS_HXX
#define POSSEL_SEQUENCESALGOS_HXX

#include "Sequences.hxx"

namespace SeqAlg
{
    /**
     * Computes the probability that randomly picking two haploytpes will be identical by descent.
     * @param seqs the sequences used to compute ehh
     */
    template<typename SeqsPtr>
    double ehh(SeqsPtr& seqs)
    {
        std::unordered_map<std::string, std::size_t> hapCounts;

        auto c = seqs->nSequences(); // c is for num of core haplotypes
        assert(c > 0);

        for(int i = 0; i < c; i++)
        {
            std::string hap = seqs->getString(i);
            if(!hapCounts.count(hap))
                hapCounts[hap] = 0;
            else
                hapCounts[hap]++;
        }

        double num{0};
        for(auto& p : hapCounts)
            num += p.second*(p.second-1)/2;
        return num/(c*(c-1)/2);
    }

    // Tajima's D-----------------------------------------

    // a1 is described in tajmad1.pdf
    double a1(std::size_t nSequences)
    {
        double sum{0};
        for(std::size_t i = 1; i < nSequences; i++)
            sum += 1.0/i;
        return sum;
    }

    double a2(std::size_t nSequences)
    {
        double sum{0};
        for(std::size_t i = 1; i < nSequences; i++)
            sum += 1.0/(i*i);
        return sum;
    }

    double b1(std::size_t nSequences)
    {
        double n = static_cast<double>(nSequences);
        return (n+1.0)/(3.0*(n-1));
    }

    double b2(std::size_t nSequences)
    {
        double n = static_cast<double>(nSequences);
        return (2*(n*n+n+3))/(9*n*(n-1));
    }

    double c1(std::size_t nSequences) { return b1(nSequences)-1.0/a1(nSequences); }

    double c2(std::size_t nSequences)
    {
        double n = static_cast<double>(nSequences);
        return b2(nSequences) - (n+2)/(a1(nSequences)*n) + a2(nSequences)/(a1(nSequences)*a1(nSequences));
    }

    double e1(std::size_t nSequences){return c1(nSequences)/a1(nSequences);}
    double e2(std::size_t nSequences){return c2(nSequences)/(a1(nSequences)*a1(nSequences)+a2(nSequences));}

    template<typename SeqsPtr>
    int tajS(SeqsPtr& s)
    {
        int segs{0};
        for(std::size_t i = 0; i < s->seqLength(); i++)
            segs += s->isSegregating(i);
        return segs;
    }

    template<typename SeqsPtr>
    double tajPI(SeqsPtr& s)
    {
        auto n = s->nSequences();
        int diff{0};
        for(std::size_t i = 0; i < n; i++)
            for(std::size_t k = i+1; k < n; k++)
                for(std::size_t allele{0}; allele < s->seqLength(); allele++)
                    diff += s->allele(i,allele) != s->allele(k,allele);

        return static_cast<double>(diff)/(n*(n-1.0)/2);
    }

    template<typename SeqsPtr>
    double tajD(SeqsPtr& seqs)
    {
        double pi = tajPI(seqs);
        double s = tajS(seqs);
        return (pi-s/a1(seqs->nSequences()))/(std::sqrt(e1(seqs->nSequences())*s+e2(seqs->nSequences())*s*(s-1)));
    }

    template<typename SeqsPtr>
    void print(SeqsPtr& seqs)
    {
        auto n = seqs->nSequences();
        std::cout << seqs->nSequences() << " sequences, " << seqs->seqLength() << " sites" << std::endl;
        std::cout << "pi: " << tajPI(seqs) << std::endl;
        std::cout << "Segregating sites: " << tajS(seqs) << "/" << seqs->seqLength() << std::endl;
        std::cout << "theta_hat[estimated from S]: " << tajS(seqs)/a1(n) << std::endl;
        std::cout << "Tajima's D: " << tajD(seqs) << std::endl;
        std::cout << "a1=" << a1(n) << " a2=" << a2(n) << " b1=" << b1(n) << " b2=" << b2(n) << std::endl;
        std::cout << "c1=" << c1(n) << " c2=" << c2(n) << " e1=" << e1(n) << " e2=" << e2(n) << std::endl;
    }
}
#endif //POSSEL_SEQUENCESALGOS_HXX
