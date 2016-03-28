//
// Created by ecorona on 3/25/16.
//

#include "gtest/gtest.h"
#include <Sequences.hxx>

TEST(GENETICMAP, GENETICMAP_DOWNLOAD_Test)
{
    for(short chr = 1; chr < 23; chr++)
    {
        Gen::GeneticMap genMap(chr);
        ASSERT_TRUE(genMap.dataExists());
    }
}

TEST(SEQUENCES, SEQUENCES_ISVALID_Test)
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

    ASSERT_TRUE(seqs->isValid());
}

TEST(SEQUENCES, setWindowByIndex)
{
    using namespace Gen;
    auto hapmap = std::make_shared<Gen::HapMapSequences>("hapmap3_r2_b36_fwd.consensus.qc.poly.chr19_ceu.unr.phased");
    hapmap->setBuild("b36");
    hapmap->setChr(19);

    hapmap->setWindowByIndex(100,101);
    ASSERT_TRUE(hapmap->seqLength() == 1);
}