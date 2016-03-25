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

    ASSERT_TRUE(seqs->isValid());
}