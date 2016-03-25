#include <iostream>
#include <vector>
#include <cmath>
#include "sel/Sequences.hxx"
#include "sel/SequencesAlgos.hxx"

using namespace std;

void test2()
{
    auto hapmap = std::make_unique<Gen::HapMapSequences>("hapmap3_r2_b36_fwd.consensus.qc.poly.chr19_ceu.unr.phased");
    hapmap->setBuild("36");
    hapmap->setChr(19);

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
    test2();
    return 0;
}