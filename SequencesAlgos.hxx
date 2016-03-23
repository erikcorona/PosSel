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
    template<typename Seqs>
    double ehh(Seqs& seqs)
    {
        std::unordered_map<std::string, std::size_t> hapCounts;

        auto c = seqs.nSequences(); // c is for num of core haplotypes
        assert(c > 0);

        for(int i = 0; i < c; i++)
        {
            std::string hap = seqs.getString(i);
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
}
#endif //POSSEL_SEQUENCESALGOS_HXX
