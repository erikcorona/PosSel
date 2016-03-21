#include <iostream>
#include <vector>
#include <cmath>
#include "Sequences.hxx"

using namespace std;

void test1(std::unique_ptr<Gen::Sequences>& seqs)
{
    seqs->addSeq("ATAATAAAAAAATAATAAAAAAATAAAAAAAATAAAAAAAA");//, "sub0");
    seqs->addSeq("AAAAAAAATAAATAATAAAAAAATAAAAAAAAAAAAAAAAA");//, "sub1");
    seqs->addSeq("AAAATAAAAATATAATAAAAAAATATAAAAAAAAAAAAAAA");//, "sub2");
    seqs->addSeq("AAAAAAAAAAAATAATAAAAAAATAAATAAATAAAAAAAAA");//, "sub3");
    seqs->addSeq("AAAATAAAAAAAATATAAAAAAATAAAAAAAAAAAAAAAAA");//, "sub4");
    seqs->addSeq("AAAATAAAAAAAAAATAAAAAAAAAAAAAAAAAAATAAAAA");//, "sub5");
    seqs->addSeq("AAAAAATAAAAATAATAAAAAAATAAAAAAAAAAAAAAAAA");//, "sub6");
    seqs->addSeq("AAAAAAAAAAAAAAATAAAAAAATAAAAAAAAAAAAAAATA");//, "sub7");
    seqs->addSeq("AAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAA");//, "sub8");
    seqs->addSeq("AAAAAAAAAAAAAAATAAAAAAATAATAAAAAAAAAAAAAA");//, "sub9");

    cout << "Checking validity: " << (seqs->isValid() ? "success" : "failed") << endl;
    seqs->print();

    constexpr std::size_t n = 100000;
    constexpr double conf = 0.05;
    auto res = seqs->sortedSamples(n, [&](auto& a, auto& b){return a->tajD() < b->tajD();});
//    for(auto& r : res)
//        std::cout << r->tajD() << std::endl;
    std::cout << conf*100 << "% [" <<
    res[std::round(  conf/2*n)]->tajD() << "," <<
    res[std::round((1-conf/2)*n)]->tajD() << "]" << std::endl;
}

int main()
{
    auto hapPtr = std::unique_ptr<Gen::Sequences>(new Gen::HaploidSequences);
    auto dipPtr = std::unique_ptr<Gen::Sequences>(new Gen::DiploidSequences);
    test1(hapPtr);
    test1(dipPtr);

    return 0;
}