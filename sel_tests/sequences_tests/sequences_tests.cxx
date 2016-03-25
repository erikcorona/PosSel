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