//
// Created by ecorona on 3/25/16.
//

#include "gtest/gtest.h"
#include <Sequences.hxx>
#include <SequencesAlgos.hxx>

TEST(SEQUENCES, tajimasD)
{
    using namespace Gen;
    using SeqPtr = std::shared_ptr<Sequence>;

    auto seqs = std::unique_ptr<Gen::Sequences>(new Gen::DiploidSequences);
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

    constexpr double err = 0.0000000001;

    EXPECT_NEAR(SeqAlg::D::tajPI(seqs), 3.8888888888888888, err); //ok
    EXPECT_EQ(SeqAlg::D::tajS (seqs),  16);
    EXPECT_NEAR(SeqAlg::D::tajD(seqs), -1.4461720765579189, err);
    EXPECT_NEAR(SeqAlg::D::a1(seqs->nSequences()), 2.8289682539682536, err);
    EXPECT_NEAR(SeqAlg::D::a2(seqs->nSequences()), 1.539767731166540, err);
    EXPECT_NEAR(SeqAlg::D::b1(seqs->nSequences()), 0.40740740740740738, err);
    EXPECT_NEAR(SeqAlg::D::b2(seqs->nSequences()), 0.27901234567901234, err);
    EXPECT_NEAR(SeqAlg::D::c1(seqs->nSequences()), 0.053921645028392084, err);
    EXPECT_NEAR(SeqAlg::D::c2(seqs->nSequences()), 0.047226772001328048, err);
    EXPECT_NEAR(SeqAlg::D::e1(seqs->nSequences()), 0.019060533801591818, err);
    EXPECT_NEAR(SeqAlg::D::e2(seqs->nSequences()), 0.0049489277698963303, err);
}