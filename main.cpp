#include <iostream>
#include <vector>
#include <cmath>
#include "Sequences.hxx"
#include "SequencesAlgos.hxx"

using namespace std;

template<typename SeqsPtr>
void test1(SeqsPtr& seqs)
{
    using namespace Gen;
    using SeqPtr = std::shared_ptr<Sequence>;

    seqs->addSeq(SeqPtr(new StrSequence("ATAATAAAAAAATAATAAAAAAATAAAAAAAATAAAAAAAA")));//, "sub0");
    seqs->addSeq(SeqPtr(new StrSequence("AAAAAAAATAAATAATAAAAAAATAAAAAAAAAAAAAAAAA")));//, "sub1");
    seqs->addSeq(SeqPtr(new StrSequence("AAAATAAAAATATAATAAAAAAATATAAAAAAAAAAAAAAA")));//, "sub2");
    seqs->addSeq(SeqPtr(new StrSequence("AAAAAAAAAAAATAATAAAAAAATAAATAAATAAAAAAAAA")));//, "sub3");
    seqs->addSeq(SeqPtr(new StrSequence("AAAATAAAAAAAATATAAAAAAATAAAAAAAAAAAAAAAAA")));//, "sub4");
    seqs->addSeq(SeqPtr(new StrSequence("AAAATAAAAAAAAAATAAAAAAAAAAAAAAAAAAATAAAAA")));//, "sub5");
    seqs->addSeq(SeqPtr(new StrSequence("AAAAAATAAAAATAATAAAAAAATAAAAAAAAAAAAAAAAA")));//, "sub6");
    seqs->addSeq(SeqPtr(new StrSequence("AAAAAAAAAAAAAAATAAAAAAATAAAAAAAAAAAAAAATA")));//, "sub7");
    seqs->addSeq(SeqPtr(new StrSequence("AAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAA")));//, "sub8");
    seqs->addSeq(SeqPtr(new StrSequence("AAAAAAAAAAAAAAATAAAAAAATAATAAAAAAAAAAAAAA")));//, "sub9");

    cout << "Checking validity: " << (seqs->isValid() ? "success" : "failed") << endl;
    SeqAlg::D::print(seqs);
}

void test2()
{
    auto hapmap = std::make_unique<Gen::HapMapSequences>("hapmap3_r2_b36_fwd.consensus.qc.poly.chr19_ceu.unr.phased");
    hapmap->setBuild("36");
    int cnt{0};
    constexpr std::size_t step = 10;
    for(std::size_t i = 0; i < hapmap->trueSeqLength() - step; i+= step)
    {
        hapmap->setWindowByIndex(static_cast<std::size_t>(0), i + step);
        double d = SeqAlg::D::tajD(hapmap);
        std::cout << d << "\t" << SeqAlg::ehh(hapmap) << std::endl;
        if(d*d > 4)
            cnt++;
    }

    std::cout << cnt << std::endl;
    std::cout.flush();
}

void test3()
{

}

int main()
{
    auto hapPtr = std::unique_ptr<Gen::Sequences>(new Gen::HaploidSequences);
    auto dipPtr = std::unique_ptr<Gen::Sequences>(new Gen::DiploidSequences);
    test1(hapPtr);
    test1(dipPtr);
//    test2();

    return 0;
}