#include <iostream>
#include <vector>
#include <cmath>
#include "Sequences.hxx"

using namespace std;

void test1(std::unique_ptr<Gen::Sequences>& seqs)
{
    seqs->addSeq(std::shared_ptr<Gen::Sequence>(new Gen::StrSequence("ATAATAAAAAAATAATAAAAAAATAAAAAAAATAAAAAAAA")));//, "sub0");
    seqs->addSeq(std::shared_ptr<Gen::Sequence>(new Gen::StrSequence("AAAAAAAATAAATAATAAAAAAATAAAAAAAAAAAAAAAAA")));//, "sub1");
    seqs->addSeq(std::shared_ptr<Gen::Sequence>(new Gen::StrSequence("AAAATAAAAATATAATAAAAAAATATAAAAAAAAAAAAAAA")));//, "sub2");
    seqs->addSeq(std::shared_ptr<Gen::Sequence>(new Gen::StrSequence("AAAAAAAAAAAATAATAAAAAAATAAATAAATAAAAAAAAA")));//, "sub3");
    seqs->addSeq(std::shared_ptr<Gen::Sequence>(new Gen::StrSequence("AAAATAAAAAAAATATAAAAAAATAAAAAAAAAAAAAAAAA")));//, "sub4");
    seqs->addSeq(std::shared_ptr<Gen::Sequence>(new Gen::StrSequence("AAAATAAAAAAAAAATAAAAAAAAAAAAAAAAAAATAAAAA")));//, "sub5");
    seqs->addSeq(std::shared_ptr<Gen::Sequence>(new Gen::StrSequence("AAAAAATAAAAATAATAAAAAAATAAAAAAAAAAAAAAAAA")));//, "sub6");
    seqs->addSeq(std::shared_ptr<Gen::Sequence>(new Gen::StrSequence("AAAAAAAAAAAAAAATAAAAAAATAAAAAAAAAAAAAAATA")));//, "sub7");
    seqs->addSeq(std::shared_ptr<Gen::Sequence>(new Gen::StrSequence("AAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAA")));//, "sub8");
    seqs->addSeq(std::shared_ptr<Gen::Sequence>(new Gen::StrSequence("AAAAAAAAAAAAAAATAAAAAAATAATAAAAAAAAAAAAAA")));//, "sub9");

    cout << "Checking validity: " << (seqs->isValid() ? "success" : "failed") << endl;
    seqs->print();

    constexpr std::size_t n = 1000;
    constexpr double conf = 0.02;
    auto res = seqs->sortedSamples(n, [&](auto& a, auto& b){return a->tajD() < b->tajD();});
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