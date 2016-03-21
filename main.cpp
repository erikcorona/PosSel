#include <iostream>
#include <vector>
#include <cmath>
#include <boost/random.hpp>

using namespace std;

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
    public:

        void addSeq(std::string seq)
        {
            seqs.push_back(Sequence(seq));
        }

        bool isValid()
        {
            for(auto& sequence : seqs)
                if(seqs[0].size() != sequence.size())
                    return false;
            return true;
        }

        // a1 is described in tajmad1.pdf
        double a1()
        {
            double sum{0};
            for(int i = 1; i < seqs.size(); i++)
                sum += 1.0/i;
            return sum;
        }

        double a2()
        {
            double sum{0};
            for(int i = 1; i < seqs.size(); i++)
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

        double c1()
        {
            return b1()-1.0/a1();
        }

        double c2()
        {
            double n = static_cast<double>(seqs.size());
            return b2() - (n+2)/(a1()*n) + a2()/(a1()*a1());
        }

        double e1()
        {
            return c1()/a1();
        }

        double e2()
        {
            return c2()/(a1()*a1()+a2());
        }

        std::size_t seqLength(){return seqs[0].size();}

        double tajPI()
        {
            int diff{0};
            for(int i = 0; i < seqs.size(); i++)
                for(int k = i+1; k < seqs.size(); k++)
                    for(int allele{0}; allele < seqLength(); allele++)
                        diff += seqs[i][allele] != seqs[k][allele];

            return static_cast<double>(diff)/(seqs.size()*(seqs.size()-1.0)/2);
        }

        bool isSegregating(int i)
        {
            char ref = seqs[0][i];
            for(auto& seq : seqs)
            {
                if (seq[i] != ref)
                    return true;
            }
            return false;
        }

        int tajS()
        {
            int segs{0};
            for(int i = 0; i < seqLength(); i++)
                segs += isSegregating(i);
            return segs;
        }

        double tajD()
        {
            double pi = tajPI();
            double s = tajS();
            return (pi-s/a1())/(std::sqrt(e1()*s+e2()*s*(s-1)));
        }

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

        Sequences sample()
        {
            static boost::random::mt19937 rng(1);         // produces randomness out of thin air
            boost::random::uniform_int_distribution<> six(0,static_cast<int>(seqs.size())-1);

            Sequences sequences;

            for(int i = 0; i < seqs.size(); i++)
                sequences.addSeq(seqs[six(rng)]);
            return sequences;
        }

        template<typename Fxn>
        std::vector<Sequences> sortedSamples(std::size_t n, Fxn sorter)
        {
            std::vector<Sequences> vec(n);

            for(int i = 0; i < n; i++)
                vec[i] = sample();
            std::sort(vec.begin(), vec.end(), sorter);
            return vec;
        }

    private:
        std::vector<Sequence> seqs;
    };


}



int main()
{
    Gen::Sequences seqs;
//               01001010101011010000000101110001100100010
    seqs.addSeq("ATAATAAAAAAATAATAAAAAAATAAAAAAAATAAAAAAAA"); // subj0
    seqs.addSeq("AAAAAAAATAAATAATAAAAAAATAAAAAAAAAAAAAAAAA"); // subj1
    seqs.addSeq("AAAATAAAAATATAATAAAAAAATATAAAAAAAAAAAAAAA"); // subj2
    seqs.addSeq("AAAAAAAAAAAATAATAAAAAAATAAATAAATAAAAAAAAA"); // subj3
    seqs.addSeq("AAAATAAAAAAAATATAAAAAAATAAAAAAAAAAAAAAAAA"); // subj4
    seqs.addSeq("AAAATAAAAAAAAAATAAAAAAAAAAAAAAAAAAATAAAAA"); // subj5
    seqs.addSeq("AAAAAATAAAAATAATAAAAAAATAAAAAAAAAAAAAAAAA"); // subj6
    seqs.addSeq("AAAAAAAAAAAAAAATAAAAAAATAAAAAAAAAAAAAAATA"); // subj7
    seqs.addSeq("AAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAA"); // subj8
    seqs.addSeq("AAAAAAAAAAAAAAATAAAAAAATAATAAAAAAAAAAAAAA"); // subj9

    cout << "Checking validity: " << (seqs.isValid() ? "success" : "failed") << endl;
    seqs.print();

    constexpr std::size_t n = 10000;
    constexpr double conf = 0.05;
    auto res = seqs.sortedSamples(n, [&](Gen::Sequences& a, Gen::Sequences& b){return a.tajD() < b.tajD();});
    std::cout << conf*100 << "% [" <<
    res[std::round(  conf/2*n)].tajD() << "," <<
    res[std::round((1-conf/2)*n)].tajD() << "]" << std::endl;
//    res.size() << std::endl;
//    cout << "Computing D: " << seqs.tajD() << endl;
    return 0;
}