//
// Created by ecorona on 3/25/16.
//

#include "gtest/gtest.h"
#include <Sequences.hxx>
#include <SequencesAlgos.hxx>

TEST(SEQUENCES, integrateEHH)
{
    auto hapmap = std::make_shared<Gen::HapMapSequences>("hapmap3_r2_b36_fwd.consensus.qc.poly.chr6_ceu.unr.phased");

    hapmap->setBuild("36");
    hapmap->setChr(6);
    Gen::GeneticMap genMap(hapmap->getChr());

    int cnt{0};
    constexpr std::size_t step = 15;
    for(std::size_t i = 0; i < hapmap->trueSeqLength() - step; i+= step)
    {
        hapmap->setWindowByIndex(static_cast<std::size_t>(i), i + step);
        std::cout << genMap.getcM(hapmap->getPos(i)) << "\t" << SeqAlg::D::tajD(hapmap) << std::endl;
    }

    std::cout << cnt << std::endl;
    std::cout.flush();
}

TEST(SEQUENCES, EHH)
{
    using namespace Gen;
    auto hapmap = std::make_shared<Gen::HapMapSequences>("hapmap3_r2_b36_fwd.consensus.qc.poly.chr19_ceu.unr.phased");
    hapmap->setBuild("b36");
    hapmap->setChr(19);
    ASSERT_TRUE(hapmap->nSequences() > 0);
    hapmap->setWindowByIndex(100,101);
    EXPECT_NEAR(SeqAlg::ehh(hapmap), 0.51336898395, 0.00000001);
}

TEST(SEQUENCES, tajimasD)
{
    using namespace Gen;
    using SeqPtr = std::shared_ptr<Sequence>;

    auto seqs = std::shared_ptr<Gen::Sequences>(new Gen::DiploidSequences);
    seqs->addSeq(SeqPtr(new Sequence("ATAATAAAAAAATAATAAAAAAATAAAAAAAATAAAAAAAA")));//, "sub0");
    seqs->addSeq(SeqPtr(new Sequence("AAAAAAAATAAATAATAAAAAAATAAAAAAAAAAAAAAAAA")));//, "sub1");
    seqs->addSeq(SeqPtr(new Sequence("AAAATAAAAATATAATAAAAAAATATAAAAAAAAAAAAAAA")));//, "sub2");
    seqs->addSeq(SeqPtr(new Sequence("AAAAAAAAAAAATAATAAAAAAATAAATAAATAAAAAAAAA")));//, "sub3");
    seqs->addSeq(SeqPtr(new Sequence("AAAATAAAAAAAATATAAAAAAATAAAAAAAAAAAAAAAAA")));//, "sub4");
    seqs->addSeq(SeqPtr(new Sequence("AAAATAAAAAAAAAATAAAAAAAAAAAAAAAAAAATAAAAA")));//, "sub5");
    seqs->addSeq(SeqPtr(new Sequence("AAAAAATAAAAATAATAAAAAAATAAAAAAAAAAAAAAAAA")));//, "sub6");
    seqs->addSeq(SeqPtr(new Sequence("AAAAAAAAAAAAAAATAAAAAAATAAAAAAAAAAAAAAATA")));//, "sub7");
    seqs->addSeq(SeqPtr(new Sequence("AAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAA")));//, "sub8");
    seqs->addSeq(SeqPtr(new Sequence("AAAAAAAAAAAAAAATAAAAAAATAATAAAAAAAAAAAAAA")));//, "sub9");

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

TEST(SEQUENCES, initialization)
{
    using namespace Gen;
    using SeqPtr = std::shared_ptr<Sequence>;

    auto seqs = std::shared_ptr<Gen::Sequences>(new Gen::DiploidSequences);
    seqs->addSeq(SeqPtr(new Sequence("ATAATAAAAAAATAATAAAAAAATAAAAAAAATAAAAAAAA")));//, "sub0");
    seqs->addSeq(SeqPtr(new Sequence("AAAAAAAATAAATAATAAAAAAATAAAAAAAAAAAAAAAAA")));//, "sub1");
    seqs->addSeq(SeqPtr(new Sequence("AAAATAAAAATATAATAAAAAAATATAAAAAAAAAAAAAAA")));//, "sub2");
    seqs->addSeq(SeqPtr(new Sequence("AAAAAAAAAAAATAATAAAAAAATAAATAAATAAAAAAAAA")));//, "sub3");
    seqs->addSeq(SeqPtr(new Sequence("AAAATAAAAAAAATATAAAAAAATAAAAAAAAAAAAAAAAA")));//, "sub4");
    seqs->addSeq(SeqPtr(new Sequence("AAAATAAAAAAAAAATAAAAAAAAAAAAAAAAAAATAAAAA")));//, "sub5");
    seqs->addSeq(SeqPtr(new Sequence("AAAAAATAAAAATAATAAAAAAATAAAAAAAAAAAAAAAAA")));//, "sub6");
    seqs->addSeq(SeqPtr(new Sequence("AAAAAAAAAAAAAAATAAAAAAATAAAAAAAAAAAAAAATA")));//, "sub7");
    seqs->addSeq(SeqPtr(new Sequence("AAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAA")));//, "sub8");
    seqs->addSeq(SeqPtr(new Sequence("AAAAAAAAAAAAAAATAAAAAAATAATAAAAAAAAAAAAAA")));//, "sub9");

    ASSERT_TRUE(seqs->allele(0,0) == 'A');
}